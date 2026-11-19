# WP-232l — two-phase particles: the interface (232l-C, round 1: the freeze candidate)

Status: **proposed** by 232l-C, 2026-11-18; binding once the orchestrator says "frozen". Authority:
ADR 232 (§4, §5.2, §5.3, §6, §8 P6, §9 phase 6, amendment "232l kickoff"), ADR-71 §3 / §6, ADR-73,
and `wp-briefs.md` "232l" (U1–U12, T0–T12). Anchors are on `origin/wp/232l-mpm-two-phase` @
`b751dcbfe`; line numbers drift. No `tests/*_reference.py` was read. Decisions are numbered **P*n***.
Section 1 is the formulation 232l-V's oracle is written from; it is self-contained.

## 1. The formulation (both routes) — symbols, signs, state, step order

**Signs.** Tension-positive Cauchy stress, compression-positive pore pressure p (excess and total
are the same unknown); σ = σ′ − α m p, m = [1,1,(1),0,0,(0)] in the Voigt order of
`tests/mpm_interface.md` §1 (LadrunoUP's split, `LadrunoUP.cpp:2235-2236`). A particle's `sig` is
**σ′** (what its material returns); total stress is never stored. Darcy flux q = −k̄(∇p − ρ_f b_f).

**Constants per two-phase body (P1).** K_f, n, ρ_f, k̄ (diagonal, = k_hydraulic/γ_w), α (default 1),
K_s (default −1: incompressible grains), stabilization (default `auto 0.25`). 1/Q̄ = n/K_f + (α − n)/K_s
(the K_s term dropped when K_s ≤ 0) — `oneOverQbar`, `LadrunoUPKernel.h:189-194`. The body's `-rho`
is the **saturated mixture** density (LadrunoUP's pinned convention, ADR-71 §3.5). Fluid drive
b_f = λ·(the `-fluidBody` vector, default: `mpm gravity`), λ the host series. Porosity is constant
(U3); a J-update would carry n per particle, n₁ = 1 − (1 − n₀)/J, and rebuild 1/Q̄ (and k̄ under a
Kozeny–Carman option) each step — not in phase 6.

**Particle state (P2).** Every particle gains `p` (committed pore pressure); mode A also `pd` (ṗ)
and `pdd` (p̈). Single-phase particles carry 0. **Objection to U4 (amendment A1):** U4 names p and ṗ;
Newmark's re-seed reads the committed accel of every DOF too (`Newmark.cpp:396-416`), and the
p-slot carries one, so p̈ must ride the particle or every step restarts Newmark's p history at p̈ = 0
(a different one-step method than LadrunoUP's).

**Shape weights.** Per particle: S_a, ∇S_a at the step-start position (the existing supports:
GIMP / linear / B-spline, renormalised in mode A, folded on B-spline faces — §3). The p weights
S^p_a, ∇S^p_a equal S_a, ∇S_a except on folded B-spline faces (§3, P9). V_p is the particle's
step-start volume: V₀ under `mpm kinematics small`, J_n V₀ under `hypo` (P6).

**Particle blocks** — the kernel helpers verbatim, (N, ∇N, dv) = (S^p or S, ∇S^p or ∇S, V_p):
`addQ(∇S, S^p, nen, nen, nd, α, V_p)` → Q_(a,i),b = α V_p ∂_iS_a S^p_b (`:87-96`);
`addH(∇S^p, k̄, V_p)` (`:104-113`); `addS(S^p, 1/Q̄, V_p)` (`:120-127`);
`addHtilde(∇S^p, 1.0, V_p)` → H̃ (`:136-145`); `addFseep(∇S^p, k̄, ρ_f, b_f, V_p)` (`:153-161`).
α_st = 0 (`off`), the raw value (manual), or `stabAlphaAuto(h, K_d, G_d, α₀)` = α₀h²/(K_d + 4G_d/3)
(`:175-178`) with **h = the grid spacing on every basis** and K_d, G_d from the particle material's
initial tangent exactly as `LadrunoUP::stabAlphaValue` reads them (`LadrunoUP.cpp:749-772`:
G = D₀(2,2) [3D: D₀(3,3)], λ = D₀(0,1), K = λ + 2G/3). S* = S + α_st H̃.

**Route E (mode A) element, per support node a, DOFs [u_a (nd), p_a] (P3):** tangent [K_u, −Q; 0, H],
damping [C_Ray, 0; Qᵀ, S*], mass [M_u, 0; 0, 0]; residual
R = [f′ − Q p − f_body ; H p − f_seep] + M a + Damp·v, with p, ṗ, p̈ the trial disp / vel / accel of
the p slot — `LadrunoUP.cpp:1360-1408`, `:1450-1503`, `:1505-1521`, `:1523-1674`, `:1682-1708`
placed one for one. f′_a = V σ′·g_a, K_u, M_u = S_a m_p I, f_body = S_a λ(m_p b + f_p) are mode A's
shipped terms (interface `mpm_implicit_interface.md` D1.6–D1.7); C_Ray = αM M_u + βK K_u from the
**solid blocks only** (ADR-71 §3.2 obligation 2). `-dynSeepage` stays off (refused if `on`).

**Route S (mode F) step** — explicit skeleton, implicit pressure (P11). With v̂_a = P_a/m_a, the
projected p̂_a (P5), lumped M⁻¹ (1/m_a on a free velocity component, 0 on a fixed or prescribed one):
1. P2G as today (`LadrunoMPMCore.h:797-871`): m_a, P_a, f_ext, f′ (effective stress only), folds.
2. Predictor (today's step 4, `:873-901`, with local damping on f′ + f_ext only):
   v*_a = v̂_a + Δt (f′_a + f_ext,a)/m_a, then the velocity conditions.
3. Pressure, on the p-active, non-drained nodes (P4, P8), drained p = p₀ moved to the right side:
   **[S* + Δt H + Δt² Qᵀ M⁻¹ Q] p = S* p̂ − Δt Qᵀ v* + Δt f_seep.** SPD. Jacobi-PCG, rel. residual
   1e-13, maxit 4N + 1000, warm start p̂ (the overlay's algorithm, `LadrunoPorousOverlay.cpp:887-982`,
   tolerance tightened from its 1e-10 for T0); non-convergence fails the step (−1, nothing changed).
4. Correction: v₁ = v* + Δt M⁻¹ Q p on free components; a_grid = (v₁ − v̂)/Δt.
5. MUSL unchanged (`:903-1013`): particle v, x from v₁; second projection; Δε from it; σ′.
6. G2P of p (P5); commit; `Domain::commit` as today.

Eliminating v₁ from M(v₁ − v̂)/Δt = f′ + f_ext + Q p and S*(p − p̂)/Δt + Qᵀv₁ + H p = f_seep
gives step 3 exactly: the scheme is IMEX (σ′ explicit, the p–u coupling and p implicit), with no
splitting error beyond that.

**Step order, mode A (route E):** today's `prepareStep` (`LadrunoMPMImplicit.cpp:704-1072`) plus:
P2G of p, ṗ, p̈ (P5) beside v, a; every grid node created with ndf = nd + 1 (P4); committed p-slot
disp / vel / accel = p̂, ṗ̂, p̈̂ (p₀, 0, 0 on a `pfix` node); stray-p and plate SPs (P4, P10); the
unconditional `domainChange()`. Newmark's `domainChanged()` re-reads all three (W1). G2P (batched,
D1.9) adds P5's update of p, ṗ, p̈.

## 2. U2 — the grid pressure DOF in mode A

- **P4 (U2; amended (defect fix 3), §7b).** In a model with any two-phase body, **every** grid node has ndf = nd + 1, p last
  (`dofMap`, `LadrunoUPKernel.h:203-216`); a model with none keeps ndf = nd bit for bit (T10). The
  node constructor's literal ndf (`LadrunoMPMImplicit.cpp:885-886`) takes nd + 1, and
  `LadrunoMPMPoint::getNumDOF` (`LadrunoMPMPoint.h:188`) returns nen·(nd + 1) for **every** particle,
  single-phase ones with zero p rows: `FE_Element::setID` fills the element's ID from each node's
  whole DOF_Group ID and refuses a count mismatch (`FE_Element.cpp:352-358`), so a single-phase
  element on an nd + 1 node must declare its p slot. **Stray p:** a node's p is *p-active* iff some
  two-phase particle's kept support has S^p_a > 0 there; otherwise the host fixes it (SP p = 0). A
  stray p couples to nothing (no single-phase particle assembles Q), so fixing it changes no answer
  and no p row reaches the solver singular. Mixed ndf was rejected: a node changing class would need
  a retire-and-recreate with a new tag, and every single-phase element would still need per-node ndf.
- **Numberer and handler** see ordinary DOFs; a fixed p is removed by Plain / Transformation, held by
  Penalty / Lagrange. **`renormalize`** (`:764-788`) drops nodes on total mass and rescales S, ∇S;
  S^p is rescaled by the same sums (B9 rule), so a dropped node takes its p with it.
- **Cost:** mode A element and memory grow by (nd + 1)²/nd² (2D 2.25×, 3D 1.78×) in a two-phase epoch.

## 3. U7, U8 — pressure conditions, Terzaghi's load, Mandel's plate

- **P8 `mpm pfix <face> <p₀>`** — p = p₀ on that grid face (p₀ = 0: drained). GIMP / linear: host SPs
  on the face nodes and their ghost layer (as `fix`, `:800-810`; the ghost layer is never reached by
  a body flush with the face: a GIMP support ends at x_p + l_p + h, open). Mode F: p eliminated there.
  On B-spline faces p₀ ≠ 0 is refused in phase 6. Impervious is natural. Drained lines inside the grid
  are not offered: every gated case drains on a grid face (U7).
- **P9 B-spline images per component (amended (defect fix 2), §7b).** Today a `fix` face folds every component of its ghosts
  (`foldedSupport`, `LadrunoMPMCore.h:462-527`; signs `:492`). New rule: a face folds **u** iff it
  carries a u condition (`fix` or `plate`), exactly as today; it folds **p** iff it folds u or has a
  `pfix` — odd under `pfix`, even otherwise (impervious). A `pfix`-only face (Mandel's drained,
  traction-free side) therefore keeps its ghost nodes as ordinary u nodes and folds only p onto the
  mirror; the ghost's own p is then stray (P4). `foldedSupport` gains a per-component target map
  (today one target per axis serves all components). Mode F folds p only inside the route-S assembly.
- **P10 `mpm plate <face> <F>` — a rigid frictionless plate, the load for both Terzaghi and Mandel.**
  F: total normal force on the plate (per unit thickness in 2D), positive into the body, × λ. The
  normal displacement component of every face-row node is replaced by one plate unknown U
  (mode A: a host-created plate node, ndf nd + 1, its other components and p fixed by host SPs, the
  load a host NodalLoad; mode F: one extra grid DOF with mass Σ folded m_a); the face-row nodes' own
  normal components become stray and are fixed. Tangential components stay free (frictionless);
  p is impervious unless the face also has `pfix` (Terzaghi's drained top). On B-spline the ghost
  row's normal component folds as an affine odd image about U, u_g = 2U − u_mirror (weights −1 on the
  mirror, +2 on U); its mass lumps onto the mirror, as today's odd fold does. GIMP / linear: exact.
  A 1D laterally-confined column under the plate is exactly Terzaghi's uniform load; the FE
  reference ties the top nodes' u_y (`equalDOF`) and loads the retained node (232l-V).
- **P10b `mpm traction <name> <face> <lo> <hi> <t: nd>`** — T5's B4 strip: consistent nodal forces
  f_a = λ t ∫_seg N_a dΓ with the linear hat along the face (GIMP / linear; refused on B-spline in
  phase 6), = FE's q·h/2 per strip node on B4's one-cell strip.
- **Ramp for T6:** T_r ≥ 20 L/c_u, c_u = √((M_d + α²Q̄)/ρ). The Terzaghi deck
  (`test_ladruno_up_element_analytic.py:216-219`: E 1e4, ν 0.2, ρ 2, K_f 4444, n 0.4, k̄ 1e-4, L 10)
  has M_d 11111, Q̄ 11110, c_u 105.4 m/s, so T_r = 2.0 s (21 transits; Tv(2 s) = 0.011).

## 4. U4, U5, U9 — projection, integrators, refusals

- **P5 (U4; projection replaced by P5a, §7b).** Projection: volume-weighted, the row-sum-lumped L2 projection at particle quadrature,
  p̂_a = Σ_p S^p_pa V_p p_p / Σ_p S^p_pa V_p over two-phase particles (ṗ̂, p̈̂ alike, mode A). Update:
  p_p ← α_p (p_p + Σ_a S^p_pa (p_a − p̂_a)) + (1 − α_p) Σ_a S^p_pa p_a, the same blend for ṗ, p̈;
  `mpm pressureBlend <α_p>`, default 1 (incremental, FLIP-like). The reference p̂ is the
  unconditioned projection, so drainage at a `pfix` node reaches the particles beside it (232f
  review item 2's rule for v). With α_p = 1 a motionless, sealed, k̄ = 0 column keeps p bit for bit;
  α_p = 0 re-smooths every step — the U4 diffusion row reports both.
- **P6 `mpm kinematics small | hypo`** (default `hypo`, fixed per epoch, both modes): `small` is
  Δε = sym Σ_a Δu_a ⊗ ∇S_a, g_a = ∇S_a, V = V₀, R = I, no geometric stiffness, x still moves. T1
  needs it (LadrunoUP `-geom linear`: mode A has no small-strain path today,
  `LadrunoMPMPoint.h:44-49`); T0 runs under it so the oracle needs no hypo kernel. **Amendment A2.**
- **P7 refusals** (loud, at `mpm analyze` or at the command; T8 tests each; messages name the rule):
  (1) a two-phase body with a symmetric-storage SOE in mode A — the existing branch
  (`LadrunoMPMImplicit.cpp:473-494`) gains the trigger "a two-phase body ([K, −Q; 0, H] with Qᵀ in
  the damping, ADR-71 §3.2)"; route S owns its SPD solve and needs none; (2) `mode quasiStatic` —
  "cannot consolidate: a static step has no S or Qᵀ rate terms"; (3) any integrator but Newmark in
  mode A (HHT, GeneralizedAlpha refused "until a test admits them"; central difference is already
  refused by tag, `:409-415`); (4) PDMY, PIMY, FluidSolidPorous — already refused by `admitMaterial`
  (`LadrunoMPMPattern.cpp:428-438`); the message gains "porous materials are not admitted (U1)";
  (5) `mpm antiLocking fbar` with a two-phase body (no tested composition; LadrunoUP's answer is
  `-formulation bbar`); (6) a non-zero σ₀ on a two-phase particle (`mpm stress`, `-s`, Contract-1
  initial stress; Q2); (7) `-dynSeepage on`; (8) n > α, K_f ≤ 0, missing `-Kf -poro -rhoF -perm`
  (LadrunoUP's parser rules, `OPS_LadrunoUP.cpp:545-575`), unknown flags fatal; (9) `pfix` p₀ ≠ 0,
  `plate` with an image conflict, or `traction` on B-spline; (10) a restart pack < 5 with a
  two-phase body (P13). **Threads:** no new code — the particle loops are serial for every material
  and mode A prints the existing serial warning (`mpm_soil_interface.md` D5); the pressure assembly
  and CG are serial.

## 5. U6 — route S: (b), the fractional step, in u-p form

**Argument.** Route (a) as the fork measured it cannot give a drained Δt. ADR-73 P0 item 3: the
rate form is "CD-unstable at any Δt"; the explicit lane is L = 0, stable to "exactly 1.000× the
[undrained] pencil"; the fixed-point-L variant is drained-stable but carries "a Δt-independent O(1)
diffusion-rate artifact (0.67 relL2)". The cause is structural: the overlay applies −Q pⁿ as a
**frozen** force, so the undrained stiffness Q S*⁻¹ Qᵀ enters the explicit operator. Route (b) puts
p^{n+1} inside the momentum update (step 4) and eliminates v₁ (step 3); the coupled volumetric
stiffness is then implicit, the explicit operator is the drained skeleton alone, and the limit is the
drained one (Huang & Zienkiewicz 1998; Book §5.5, "Δt ≤ h/c_s", ADR-71 §3.3). **Berkeley's
structure, not its equations — amendment A3.** `MPMSemiImplicitTwoPhase` (geomechanics/mpm @
`b6f2edf69`, MIT) runs predictor → Poisson → pressure → correction force
(`mpm_semi_implicit_twophase.tcc:335-353`), but on **two velocity fields with drag and
incompressible constituents**: the Poisson right side is −D_s v_s − D_l v_l
(`assembler_eigen_semi_implicit_twophase.tcc:228-229`), drag `:33-49`. That is u-U (zone III), out
of U1, with no storage S. Route S keeps Berkeley's step structure on ADR-71's single-velocity u-p
equations. Its projection parameter (β, default 1 = incremental, `mpm_semi_implicit_twophase.h:101`,
used at `assembler…:296-317`) has no counterpart: the elimination is exact, so solving for p or
for p − p̂ gives the same p^{n+1}.
- **Δt limit.** Δt_S = 2/ω_max of the drained pencil M_L⁻¹K_d (lumped grid mass, the particle
  tangents, constrained components removed): `mpm criticalTimeStep drainedPencil` (power iteration,
  200 iterations) is T7's yardstick, the discrete pencil as ADR-73 P0 item 3 rules. The run-time
  refusal stays today's estimator, Δt ≤ s·h / max_p(c_d + |v_p|), c_d = √(M_d/ρ_mix) — no fluid term
  (`LadrunoMPMCore.h:745-760`). `mpm criticalTimeStep undrainedPencil` (D_u = D + α²Q̄ m mᵀ at each
  particle) is T7's reported undrained limit. The H term is implicit: no diffusion limit.
- **Stabilization.** S* carries the body's `-stab`, as in route E. The elimination adds
  Δt² QᵀM⁻¹Q ≈ (Δt/Δt_S)²·h²/M_d, a McGann-type α₀ ≈ C² Laplacian; whether it alone suppresses B4's
  checkerboard is reported (R: route S with `-stab off`), the gate runs the default.
- **Reuses from the overlay:** the kernel blocks, the rate-form bookkeeping of p̂, the Jacobi-PCG
  algorithm and its drained-row handling — re-written as engine code on CSR over p-active nodes;
  the overlay object is not reused (FE cells only, `.cpp:229-231`). Lives in a new
  `LadrunoMPMTwoPhase.h`, called from `advance()` between step 4 (`:873`) and step 5 (`:903`).

## 6. U11, U12 — output and restart

- **P12 (U11).** `232a_mpm_contracts.md:636`: "A new reserved field (e.g. a pore pressure in phase 6)
  is additive … A minor-2 reader that does not know it skips it rather than refusing it." The fork's
  validator checks the catalogue only for names it knows (`ladruno_format.py:1343-1348`). So:
  `PORE_PRESSURE`, COMPONENTS `p`, DIMENSION `F/L^2`, DATA_TYPE 0, compression positive; 0 for
  single-phase particles; no version bump. `STRESS` of a two-phase particle is σ′ (a format-doc line).
  Edits: `ParticleWriter::components` (`LadrunoMPMIO.cpp:50-64`), `PARTICLE_FIELDS`
  (`ladruno_format.py:1015`), §2.6 of the contracts doc. apeGmsh's reader is queued with T, not here.
  Queries: `mpm particleData p | pd | pdd`; `mpm particlePressure <p..> [<pd..> <pdd..>]` sets them
  (flat, storage order, as `particleVelocity`); `mpm info twoPhase <body>` →
  [route 1 E | 2 S, K_f, n, ρ_f, k̄.., α, K_s, stabMode, α₀ or α, α_st(h), pressureCutoff]
  (the last entry added by ruling c84f71fb3).
- **P13 (U12).** Mode-A pack **5** (`MPM_PACK_VERSION = 4`, `LadrunoMPMImplicit.cpp:1533`): per
  particle p, pd, pdd; per body the fluid constants; α_p, kinematics, pfix, plates (with the plate
  node's tag), tractions. A pack < 5 with a two-phase body is refused; without one it restores as
  today. Mode F keeps refusing a database save with particles.

## 7. T0, the debug hooks, files, cases, cost

**T0 (232l-V's oracle, from §1–§5 only; ElasticIsotropic, `kinematics small`):**
(i) `mpm debug particleBlocks <p> <path>`: arms particle p's element at the committed state, sets
its trial support values to the pattern value_j = s_k sin(1.7 j + 0.3) (j the element DOF index,
s_k per kind: u, p, v, ṗ, a, p̈), dumps nodes, coordinates, S, ∇S, K, C, M and R incl. inertia, then
restores; (ii) three mode-A Newmark steps (`algorithm Linear`, a 2×2-cell GIMP grid filled flush,
so `renormalize` drops nothing, a small sine velocity so supports move): particle x, v, a, σ′, p,
ṗ, p̈ after each step; (iii) ten route-S steps on a 4×4 grid (gravity, `pfix` top, rollers, a
non-uniform initial p): the same fields plus the grid p. 232l-V chooses the numbers and writes
them as `t0_*` blocks. C++ matches each at rel 1e-9. **T1:** `mpm debug assemble <path>` dumps the
Newmark first-iteration tangent and unbalance at the armed state, each DOF labelled (x, y[, z],
dof), nothing committed; the common state is σ′ = 0, F = I, given particle (v, a, p, ṗ, p̈), the FE
nodal state being its P5 projection (232l-V computes it); grid nodes are read after
`mpm analyze 1 dt` with `nodeDisp` / `nodeCoord`. LadrunoUP side: `-lumped`, `-stab` as the leg,
`-dynSeepage off`, linear geometry.

**Cases file `tests/mpm_twophase_cases.txt`** (232l-V creates it): the format of
`mpm_implicit_cases.txt`; new keys `twophase <body> Kf n rhoF perm.. [alpha] [Ks] [stab ..]`,
`pfix <face> <p0>`, `plate <face> <F>`, `traction ..`, `kinematics small|hypo`, `pblend <α_p>`,
`pinit <body> uniform <p> | list`; observables `p_err_max@Tv`, `U_err@Tv`, `p_center@t`, `cb_lap`,
`settle0_ratio`, `hydro_p_err`, `hydro_syy_err`, `stable_steps`. `kind reference` blocks
(232l-V: data under `tests/data/mpm/twophase/`) then `kind run | refusal | protocol` blocks (mine).

**Test files (all `tests/test_ladruno_mpm_twophase_*.py`, harness
`tests/ladruno_mpm_twophase_harness.py`):**

| File | Rows | Zone-A cost |
|---|---|---|
| `_oracle` | T0 i–iii | zone_a, ~5 s |
| `_fe_identity` | T1 Q4 one cell + 2×2, H8 one cell; the T11 red oracle | zone_a, ~6 s |
| `_consolidation` | T2, T3, T6 | zone_a: Terzaghi to Tv 0.1, gimp, both routes (~10 s); rest zone_b |
| `_hydrostatic`, `_checkerboard` | T4, T5 | zone_b |
| `_timestep` | T7 | zone_a, ~5 s (2×2000 tiny mode-F steps) |
| `_refusals`, `_restart`, `_engines` | T8, T9, T12 | zone_a, ~3 / ~5 / ~10 s (classic Tcl subprocess) |

Zone-A total ≈ 45 s. **Cost per leg** (one core, 2D; mode-A step ≈ 1–2 ms per 40 particles plus
the rebuild; route-S step ≈ 0.1 ms per 100 particles): T2 gimp (40 particles, 4500 steps of 0.02 s)
~10 s, bspline2 (80) ~20 s, linear ~10 s; T3 (≈ 400 particles, ≈ 2000 steps) 0.5–1 min per leg;
T4 (to Tv ≈ 2, ≈ 18 000 steps) ~1 min; T5 B4 (≈ 400 particles, 100 steps of 0.01 s) < 1 min per
leg × 3; T6 Terzaghi (≈ 13 400 route-S steps at s = 0.5) ~5 s, Mandel 1–5 min; T7 seconds. The
whole battery fits under 30 min; the warrant run goes to esmeralda as the brief says.

## 7a. Addendum for the oracle (232l-C, 2026-11-18; answers the orchestrator's five questions)

Nothing in P1–P13 moves. This adds the details §1 left implicit. An oracle written from §1 and
§7a alone must reproduce T0 to rel 1e-9.

**Q1. The `mpm debug particleBlocks` dump.** The command is `mpm debug particleBlocks <p> <path>
<s_u> <s_p> <s_v> <s_pd> <s_a> <s_pdd>`, mode A only, with p the 0-based storage index of
`mpm particleData`. It writes one JSON object to `<path>`:
`{"particle": p, "id", "ndm", "ndf", "nen", "lambda", "V", "m", "alpha_st", "nodes": [tag..],
"coords": [[x, y(, z)]..], "S": [nen], "dS": [nen][ndm], "Sp": [nen], "dSp": [nen][ndm],
"trial": {"u": [nen][ndm], "p": [nen], "v": [nen][ndm], "pd": [nen], "a": [nen][ndm],
"pdd": [nen]}, "K": [ndof][ndof], "C": [..], "M": [..], "R": [ndof], "R_static": [ndof]}`.
Matrices are row-major lists of rows; every number uses 17 significant digits.
- **Element DOF order.** Node-major, P4's `dofMap`: element DOF j = ndf·k + c, where k is the
  support slot, c = 0..ndm−1 are u_x, u_y(, u_z) and c = ndm is p (ndf = ndm + 1 in a two-phase
  epoch). The single scalar index j also places the pattern values below.
- **Support-node order.** k runs over the tensor product of the per-axis supports. On each axis
  the grid index rises (the scan i = c−1 … c+2, with S = dS = 0 entries dropped,
  `LadrunoMPMCore.h:259-276`), and the product runs x fastest, then y, then z
  (`particleSupport`, `:282-302`). Nodes removed by `renormalize` are absent; the survivors keep
  their order. A plate unknown (P10), when present, is the last slot. `nodes` and `coords` are in
  that order, so an oracle may also map by coordinates.
- **Pattern values.** For every element DOF j:
  `disp_j = s_kind sin(1.7 j + 0.3)`, `vel_j = s'_kind sin(1.7 j + 0.3)` and
  `accel_j = s''_kind sin(1.7 j + 0.3)`. The kind is u or p, by c. So (s_u, s_p) scale the
  displacement and p slots, (s_v, s_pd) the velocity and ṗ slots, and (s_a, s_pdd) the
  acceleration and p̈ slots. They are written to the support nodes' trial disp, vel and accel.
  The `trial` block repeats them. The previous trial state is restored afterwards and nothing is
  committed.

**Q2. Yes: K, C, M and R are all formed from the trial state after the pattern values are set.**
- **Sequence.** The element is armed at the committed particle state, exactly as a step's
  `prepareStep` does it (supports at the committed x; F_n, σ′_n, V_n). The pattern values are set,
  then `update()` runs: Δε = sym Σ_k u_k ⊗ ∇S_k under `small`, then the material's trial strain.
  Then K = `getTangentStiff`, C = `getDamp`, M = `getMass`,
  R = `getResistingForceIncInertia` = R_static + M·a + C·v, and R_static = `getResistingForce`.
- **u–u blocks under `small`.** K_u = V₀ Bᵀ D B, with D the material's `getTangent()` at the trial
  strain, B in the engineering-shear Voigt order of §1 and B built from ∇S. There is no
  geometric term. M_u = diag(S_k m_p) on every u DOF of slot k (lumped; S_k = W_k, the unsigned
  weight). C_u = αM M_u + βK K_u.
- **p blocks.** p rows and columns are §1's blocks: −Q in K, Qᵀ and S* in C, H in K, zero in M.
- **R terms.** f′_k = V₀ σ′(trial)·∇S_k. f_body,k = S_k λ (m_p b + f_p). f_seep uses
  b_f = λ·(`-fluidBody` or `mpm gravity`). λ is the host's current factor: 1 before the first
  step, otherwise the factor of the last accepted step. *(Superseded for π by the Fable audit,
  finding 3, wp-briefs "232l": π is evaluated at the drive of the step being taken; see §7b.)*
- **α_st.** It uses h = the grid's cell size and D₀ = the material's `getInitialTangent()`, by
  §1's rule.
- **Single-phase particle in a two-phase epoch.** Its p rows and columns are zero.

**Q3. GIMP weights after `renormalize` in T0 (ii).**
- **The procedure** (`LadrunoMPMImplicit.cpp:726-789`):
  1. Supports by `gimpFactor` (`LadrunoMPMCore.h:141-159`), with r = x_p − x_i and l_p = h/(2·ppc):
     - |r| < l_p: S = 1 − (r² + l_p²)/(2h l_p);
     - l_p ≤ |r| < h − l_p: S = 1 − |r|/h;
     - h − l_p ≤ |r| < h + l_p: S = (h + l_p − |r|)²/(4h l_p);
     - otherwise 0.
     The tensor product is S = S_x S_y and ∇S = (S′_x S_y, S_x S′_y).
  2. Nodal mass m_a = Σ_p S_pa m_p over every support node, ghosts included.
  3. A node is kept iff m_a > 0.05 · max_p m_p.
  4. Per particle, over its kept nodes: σ = Σ S, γ = Σ ∇S, S′ = S/σ and ∇S′ = (∇S − S′γ)/σ.
- **What renormalize does on the flush-filled 2×2 grid.** Every particle centre lies inside, and
  its support ends at |r| < h + l_p (open), so a body whose particle domains are flush with the
  faces reaches no ghost node. Every real node has m_a ≥ m_p (a corner holds
  Σ = (0.75 + 0.25)² m_p at ppc 2), so nothing is dropped. Then σ = 1 and γ = 0 to round-off, and
  renormalize is the identity to round-off.
- **Particles that cross a face.** A ghost node then gets S = δ²/(4h l_p), where δ is the
  overshoot, far below the cutoff. It is dropped and the formula above rescales the rest. The
  oracle should implement steps 1–4 literally, not assume the identity: the u stays free on every
  face without a `fix`.
- **S^p.** On GIMP S^p = S′ and ∇S^p = ∇S′.

**Q4. Route S's M⁻¹ rule.** M⁻¹ is diagonal over velocity components (a, i).
- **Free components.** M⁻¹ = 1/m_a. m_a is the lumped grid mass after the B-spline mass folds
  (ghost mass on its image, `:858-864`).
- **Zero entries.** A component gets 0 when it is:
  - fixed by `fix` (face node or its ghost layer);
  - prescribed by `gridvel`;
  - on an inactive node (m_a ≤ c_m · max m_p);
  - stray (P10: a plate face-row node's normal component).
  Its v₁ is then the conditioned predictor value. Its Q row still enters Qᵀv* on the right side,
  so a prescribed velocity drives the pressure. It never receives the correction Δt M⁻¹ Q p.
- **The plate unknown U.** U is one more velocity component with M⁻¹_U = 1/M_U, where
  M_U = Σ m_a over the face-row nodes (ghost-row mass stays on the mirrors, A4). Its Q row is
  Q_U,b = α Σ_p V_p S^p_pb ∂_n S′_{n,U}. S′_{n,U} is P10's plate weight of the normal component:
  the face-row S, plus 2× the ghost-row S on B-spline. Its predictor carries the force λF.
- **Folded B-spline images.** Slaved ghosts carry no components of their own. Their Q rows are
  folded into their images with P9's signs, since Q is assembled from the folded S′_{i,a}
  and S^p_b.
- **Walls.** `mpm wall` (frictional, conditional) is refused with a two-phase body in route S.
- **p-activity in route S** (orchestrator's ruling, 3d3699a50; it replaces "p-active" in §1
  step 3 and P4 for mode F).
  - **The cutoff.** A node is p-active iff two things hold:
    - m̃_a > c_m · m_max, where:
      - m̃_a is the folded nodal mass, × 2^k on a B-spline face layer
        (`imageCutoffMass`);
      - m_max = max over every particle, single- and two-phase, of m_p;
      - c_m is route S's own knob, `mpm pressureCutoff <c>` (default 5e-2, the value of mode
        A's test at `LadrunoMPMImplicit.cpp:741-753`; reported by `mpm info twoPhase`). It is
        never moved by `mpm massCutoff`, which keeps setting the u-side cutoffs only (mode F
        1e-10, mode A 5e-2). This is ruling c84f71fb3, and T8 checks that `massCutoff` leaves
        it unchanged. *(Scoped to route S, mode F, by the Fable audit, finding 4: in mode A
        p-active ⊆ u-kept, so a two-phase body refuses `massCutoff` > `pressureCutoff`; see §7b.)* Mode A's p rows follow `renormalize` (P4) and have no separate knob;
    - some two-phase particle has S^p > 0 there.
  - **Renormalisation.** Each two-phase particle's S^p and ∇S^p are renormalised over its
    p-active nodes by Q3's formula, using B9's rule on B-spline images.
  - **A dropped node** is absent from the pressure system, the projection and G2P. It is never
    treated as drained.
  - **The u side** keeps mode F's 1e-10 cutoff: Q's u rows use the unrenormalised ∇S′, its p
    columns the renormalised S^p.

**Q5. The first-step re-seed (mode A, step 1, from the initial particle state).**
- **Particle state before step 1.**
  - x, v, a: v from `mpm velocity` / `particleVelocity`, a from the new
    `mpm particleAcceleration <flat>` (default 0);
  - p, ṗ, p̈: from the new `mpm particlePressure <p..> [<pd..> <pdd..>]` or
    `mpm porePressure <body> uniform <p>` (default 0);
  - σ′ = 0, F = I, V = V₀.
- **The grid, built from it with the Q3 weights after renormalize.**
  - Mass: m_a = Σ W m.
  - Velocity and acceleration: v̂_a,i = Σ_p S′_{i,pa} m_p v_p,i / m_a, and â the same with a
    (on GIMP / linear, S′_i = S′).
  - Pressure: over the two-phase particles only, p̂_a = Σ S^p V p / Σ S^p V, with
    V = the step-start volume (V₀ under `small`); ṗ̂ and p̈̂ use the same weights. Nodes without
    two-phase support get 0.
- **Committed node state written before `domainChange()`.**
  - u slots: disp 0; vel v̂ and accel â on free components; 0 and 0 on a fixed component;
    v_prescribed and 0 on a `gridvel` component.
  - p slot: disp p̂, vel ṗ̂, accel p̈̂; p₀, 0, 0 on a `pfix` node; 0, 0, 0 on a stray-p node.
  - Newmark re-reads all of these (W1).
- **The elements' FLIP references** are the unconditioned v̂, â, p̂, ṗ̂ and p̈̂ (P5).
- **Every later step** repeats this from the committed particle state.
- **The T1 common state is this re-seed**, and the FE nodal state is its projection: u = 0,
  v = v̂, a = â, p = p̂, ṗ = ṗ̂, p̈ = p̈̂.

**`mpm debug assemble <path>`** (T1).
- **Requires** a Transient Newmark analysis.
- **Runs** `prepareStep` as step 1 would, then the integrator's `newStep(dt)` with the dt given as
  a third argument (`mpm debug assemble <path> <dt>`), which sets the Newmark predictor, then
  `formTangent` and `formUnbalance`.
- **Writes** `{"A": [n][n], "b": [n], "dofs": [[nodeTag, x, y(, z), dof]..]}`, where dof is
  0-based and the equation order is the SOE's.
- **Then reverts the Domain**, so nothing is committed.

**Case notes (no interface change).** T3 and T6 use a 5 s linear ramp, which replaces T3's step
load and P10's 2 s Terzaghi ramp; it is still ≥ 20 transits. T4 runs at h = H/20.

## 7b. Amendments after route E on `fd81cbdd3` (232l-C, 2026-11-18; proposed for the freeze)

The owner's decision and the ruled defect fixes are in `wp-briefs.md` "232l", block "Route E on
build fd81cbdd3". P5a replaces P5's projection. Fix 2 amends P9; fix 3 amends P4. Everything else
in P1–P13 and §7a stands.

### P5a — the particle pressure carrier: FLIP increments on a consistent L2 projection

**The cause, from the diagnostic (`13024db9f`).** P5 projects with the row-sum-lumped matrix,
p̂ = M_L⁻¹ Pᵀ V p_part, where:
- P is the particle-by-node matrix of S^p;
- V = diag(V_p);
- M_L = diag(Pᵀ V P 1).

Lumping makes projection followed by interpolation lossy: M_L⁻¹ Pᵀ V P ≠ I. Two things follow:
- **PIC (α_p = 0) is hopeless.** It applies that smoothing on every step, and drained T2 went to
  −0.83 p₀.
- **FLIP (α_p = 1) settles on the wrong particle values.** It stops only where the grid increments
  vanish, i.e. p̂ = g for the grid solution g, including g_c = p₀ at a pinned node. Its particles
  therefore satisfy M_L⁻¹ Pᵀ V p_part = g. Since M_L⁻¹ Pᵀ V P ≠ I, that solution is
  p_part = P g + r with a residual r ≠ 0, whose job is to undo the lumping. Near a drained face this
  is the measured residual: the top particle at −0.044 p₀ where the field is +0.045 p₀. It is
  70–75 % grid-visible, through the pinned nodes. The grid field itself is accurate, which is the
  diagnostic's (ii).

**The scheme.** For the p-active nodes of the step (fix 3), each of x ∈ {p, ṗ, p̈}:
1. **Projection (P2G).** Solve M x̂ = Pᵀ V x_part, with the consistent matrix
   M_ab = Σ_{two-phase p} V_p S^p_pa S^p_pb. S^p is the renormalised p weight (fix 3), with the
   B-spline folds of P9 (and fix 2's affine part, see below). This is the least-squares fit of the
   particle field by the grid field (Burgess, Sulsky & Brackbill 1992, the consistent-mass FLIP
   transfer; Wallstedt & Guilkey 2011, weighted-least-squares particle-to-grid reconstruction).
   - The solve: Jacobi-PCG to a relative residual of 1e-14, warm-started from 0. The matrix is SPD
     whenever the particles of each p-active node's support span its shape functions, which holds
     for ppc ≥ 2 on gimp, linear and B-spline.
   - On a solve failure the step is refused, with nothing changed.
   - x̂ is the unconditioned reference. The committed nodal state takes x̂, except that pinned nodes
     take p₀, 0, 0 (§7a Q5, unchanged).
2. **Update (G2P).** Unchanged from P5:
   x_p ← α_p (x_p + Σ_a S^p_pa (x_a − x̂_a)) + (1 − α_p) Σ_a S^p_pa x_a, with α_p = 1 the
   default.

**Why the residual cannot survive.** P Pᵀ-type projections are exact here: M⁻¹ Pᵀ V P = I. Split
any particle field into p_part = P g + r_vis + r_inv, where r_vis = P M⁻¹ Pᵀ V r and Pᵀ V r_inv = 0.
- **At a free node** the next projection returns g + M⁻¹ Pᵀ V r_vis. The visible part becomes grid
  history and is evolved by the u-p equations: it is no longer a residual.
- **At a pinned or prescribed node c** the grid holds p₀, so the increment there is
  p₀ − x̂_c = −(M⁻¹ Pᵀ V r)_c. FLIP subtracts exactly that component from the particles. After one
  step the residual has no component the pinned node can see.
- **The fixed point.** FLIP stops only at x̂ = g, which now means Pᵀ V (p_part − P g) = 0. The
  particles carry the grid field exactly, plus at most a V-orthogonal r_inv.
- **The invisible part.**
  - It never enters a solve, since its projection is zero by definition.
  - It is never created by an increment, because increments lie in the range of P. The one source
    is support change under motion: the new step's projection drops the part of p_part outside
    the new range. That is O(|Δx|/h) per step, a fraction of the step's pressure change.
  - It is zero whenever the initial field is zero or grid-representable, which covers every gated
    case (zero excess p).
  - It is reported, never gated: the row "residual" of T2p/T4p, now
    `p_p − Σ S^p_pa p_a` after G2P.
  - Damping it would break the sealed-column invariant below for a non-representable initial
    field, so P5a does not.
- **XPIC** (Hammerquist & Nairn 2017) and **APIC-type** transfers (Jiang et al. 2015) target noise
  in velocity's null space. They do not address this residual, which is lumping bias, not noise,
  so neither is adopted.

**The sealed-column invariant (guard) holds bit for bit, as P5 claimed.** With no motion, k̄ = 0 and
ṗ = p̈ = 0, Newmark's p predictor equals x̂ exactly and the right side of every p row is exactly
zero. The increment x_a − x̂_a is therefore exactly zero, and the particles keep p bit for bit,
whatever x̂ is.

**Singularity guard (frozen, wp-briefs "Freeze of P5a").**
- **(a) The solve refuses the step loudly** (rule named, nothing changed) in either of two cases:
  - Jacobi-PCG does not reach a relative residual of 1e-14 within 4 n + 1000 iterations;
  - any PCG curvature dᵀ M d falls to 1e-14 · max diag(M) · dᵀ d or below.

  A zero diagonal entry is refused the same way. The projection runs before any grid node is
  created or retired, so a refused step leaves the model unchanged.
- **(b) T8 gains a rank-deficient case:** bspline2 with one particle per cell in a p-active region.
  It must be refused, not solved.
- **(a) amended (ruling on `7b7cdc9cf`): the probe.** Each projection also runs one probe PCG on
  b_k = sin(1.7 k + 0.3) over the p-active rows in storage order, with the same tolerance, the same
  iteration cap and the same curvature test. A failed probe refuses the step loudly (rule named,
  nothing changed).
  - **Why it is needed.** The data's right side Pᵀ V x always lies in M's range, so PCG converges
    on a singular M without failing. The probe does not lie in that range.
  - **Answers are unchanged.** A nonsingular M passes the probe, and its solution is discarded.
  - **Cost.** One more PCG per step, about one third of the projection's three solves.
- **p-activity, implementation fix (ruling on `7b7cdc9cf`).** A node is p-active iff two things
  hold:
  - its two-phase nodal mass passes `pressureCutoff`;
  - some two-phase particle has S^p > 0 there (ruling `3d3699a50` in full).
- **Renormalisation (implementation note, disclosed for ruling).** The p renormalisation of P4
  amended treats a kept entry whose p weight is folded away (S^p = 0, e.g. a ghost behind a
  drained-only face) as retained: its W counts in σ and γ. Renormalisation therefore runs only when
  a p-reached node (S^p > 0) is dropped. Read literally, the text would rescale every particle beside
  a drained-only B-spline face by Σ W over its p-active nodes alone, which is a spurious scaling.
  GIMP and linear, where nothing folds, are identical under both readings.

**Cost.** One sparse SPD solve per step with three right-hand sides, over the p-active nodes. The
pattern is that of H.

**Both routes.** Route S (§5) projects p̂ the same way. Its step 3's S* p̂ is the consistent
projection.

**What the oracles must change (232l-V re-derives from this text).**
- **T0 (ii) and T0 (iii):** replace the lumped projection by solving M x̂ = Pᵀ V x_part for
  x ∈ {p, ṗ, p̈}. M is over the p-active nodes, which include the pinned ones. Use the
  renormalised S^p of fix 3 and the particles' step-start V_p. The committed p state and the G2P
  rule are unchanged. A dense solve is fine at these sizes.
- **T0 (i)** is unaffected: it uses the committed state before any step.
- **T1** is affected: its FE nodal p, ṗ, p̈ are "the P5 projection" of the particle values, so
  t1_fe.json is regenerated with the consistent projection.

### P9, amended (defect fix 2): the affine even image on an impervious face with a fluid drive

- **The flaw.** An impervious B-spline face whose normal component of the fluid drive is non-zero,
  ρ_f b_f·n ≠ 0, cannot carry zero flux with an even image. Zero flux needs ∂p/∂n = ρ_f b_f·n
  (Darcy, q = −k̄(∇p − ρ_f b_f)), and the even image imposes ∂p/∂n = 0.
- **The image.** With d the face's axis, the ghost coefficient is
  c_g = c_m + ρ_f λ b_f,d (X_g,d − X_m,d). This is the image of p = p_even + ρ_f λ b_f,d (x_d − x_face),
  and quadratic and cubic B-splines reproduce the linear part exactly.
- **Where the known term goes.** The interpolated pressure gains
  π_p = Σ_{g folded} S_g(x_p) ρ_f λ b_f,d (X_g,d − X_m,d), and ∇π_p likewise.
  - In the u rows, the coupling force gains α V_p ∇S_a π_p.
  - In the p rows, H's flux gains V_p (∇S^p_b · k̄ ∇π_p).
  - The projection fits x_part − π_p, and G2P adds back π at the step's end. *(Superseded by the q + π invariant below: G2P alone does not add
    the change of π with particle motion.)*
  - π̇ (from λ̇) is neglected; it is exact for a constant drive, as in T4.
- **Patch check (before T4p):** a bspline2 hydrostatic column with k̄ = 1, run 200 steps of dt 1 s.
  The base particle's p equals γ_w(H − y) to 1e-10 γ_w H.
- **A′ (frozen, ruling `f69e450dd`): π's invisible part, compensated once.**
  - **The gap.** The field is p = Σ S^p c + π, and π is not in the range of P. A particle field that
    starts at zero (or at any value) therefore carries the invisible remainder −π_⊥. P5a's
    projection never removes it: on build `d06695367` the grid coefficients were exact to 2e-15,
    but the patch check's base particle missed by 1.0e-3 γ_w H.
  - **The step.** At the first projection at which π ≠ 0, each two-phase particle's p gains only
    π's own invisible part, π_⊥ = π − P M⁻¹ Pᵀ V π, on that projection's p-active set, M and S^p.
    - It is applied after every solve of that projection has succeeded. If the step is then
      refused, nothing is applied.
    - Because Pᵀ V π_⊥ = 0, it leaves that projection's x̂ unchanged.
    - ṗ and p̈ are unchanged.
  - **Never again.** G2P adds the change of π in full, so the residual stays constant. *(Annotation, build
    `89718fdad`: it did not. G2P added π_end − π_start along the element's own fold, but the
    change of π with particle motion was not carried into q; the patch base missed by 1.2e-3. The
    q + π invariant below replaces this bullet's mechanism; A′ itself stands as its Init.)* A drive that
    ramps up from 0 never triggers it, because π = 0 at the first projection.
  - **State.** The flag (applied or not) and max |π_⊥| are in the restart pack (version 5), so a
    restart never re-applies it. `mpm reset` clears both with the particles.
    `mpm info twoPhase` reports max |π_⊥| applied as its last entry, after `pressureCutoff`.
  - **What it leaves alone.** GIMP, linear and every case without a normal fluid drive on an
    impervious B-spline face have π = 0, and are untouched. That includes T0, T1 and the sealed
    column. User content in the particles' representable part is never erased (A, which would
    have, was declined).
- **The "field at the particle points" readout** (the T2p/T4p reported column) includes π.
- **The q + π invariant (frozen, ruling on `89718fdad`).** One invariant; the π_c bookkeeping is
  accepted as one implementation of it.
  - **Invariant.** Each two-phase particle carries q_p, the coefficient-space part, and
    p_p ≡ q_p + π(x_p, drive), evaluated at the CURRENT position and drive wherever p is used:
    output, queries, and the projection's right side.
  - **q is FLIP-carried** (P5a in plain space): c = M⁻¹ Pᵀ V q, and G2P adds Σ S^p (c_new − ĉ).
    The PIC part (α_p < 1) is Σ S^p c_new in q-space.
  - **Init.** q = p_init − (P M⁻¹ Pᵀ V π(0))_p. This is A′.
  - **Implementation (π_c).** The particle stores p and π_c, the π of its last evaluation;
    q = p − π_c at every projection, motion included.
    - The projection's right side is V (p − π_c); at the epoch's first projection it is
      V (p_init − π), and A′ then makes q = p_init + π_⊥ − π.
    - The element's G2P works in q-space and returns p − π_c + π_c (the PIC part adds π_c back).
    - After every accepted step the host re-evaluates π at the new positions and drive (the same
      supports, p-activity and S^p as a projection, no side effect on the model) and sets
      p ← p + π − π_c, π_c ← π. Each projection does the same before it fits (a no-op unless the
      drive or the positions changed since).
    - π_c joins the restart pack (version 5, after max |π_⊥|); `mpm reset` clears it.
  - **Rates.** ṗ and p̈ still neglect π's rate.
  - **Where π = 0** (GIMP, linear, no normal fluid drive on an impervious B-spline face) nothing
    changes: π_c = 0 and q = p.
  - **Self-check (harness, numpy-free).** On the patch column, after every step,
    max_p |(p_p − π(x_p)) − (P c)_p − r_user,p| ≤ 1e-12 γ_w H, with r_user = 0 here; c is the grid
    p read after the step, π and S^p from the particle's block dump. *(Ruled at `0fc5a3f13`: the
    self-check is the bookkeeping identity |p_p − π(x_p, drive) − q_p| ≤ 1e-12, with q's FLIP
    carry and π_c = π(x_p) checked; q − P c is the motion residual, reported.)*
- **Implementation notes after the Fable carrier audit (wp-briefs "232l"; frozen text unchanged).**
  - **Renormalisation (ruling `bb1421559`).** The u side's mass-cutoff renormalisation rescales
    the u components only. The p weights and fix 2's coefficients follow fix 3 on the raw
    weights, and a folded-away entry is retained even when the u side dropped it.
  - **Finding 1.** On a face with a u fold and `pfix`, the ghost's W folds onto the mirror only
    (a dangling `else` once added it to the ghost too: ΣW = 1.03125).
  - **Finding 2.** A′'s π_⊥ and the new π_c are computed at the projection but committed (p, π_c,
    the applied flag, max |π_⊥|) only when the step is accepted: G2P leaves q_new + π_c, and the
    accepted step adds α_p π_⊥ (first step only) + π(x_new) − π_c. A refused step and the debug
    hooks change none of them. Setting p by command (`particlePressure`, `porePressure`) clears
    the flag, so A′ re-arms for the new field.
  - **Finding 3.** π is evaluated at the drive of the step being taken (the series value that
    step's `applyLoad` uses), in the projection's step and after it. A′'s π(0) is π at the epoch's
    initial drive, λ(t₀): a drive that ramps from 0 needs no A′ (max |π_⊥| = 0), and the first
    step's π_c is that π(0).
  - **Finding 4.** In mode A, a two-phase body refuses `massCutoff` > `pressureCutoff`, naming the
    rule. With massCutoff ≤ pressureCutoff, p-active ⊆ u-kept holds by itself.

### P4, amended (defect fix 3): mode A's p-activity is route S's rule

- **The rule.** A node's p is active only if its two-phase nodal mass exceeds `pressureCutoff`
  × m_max:
  - the two-phase nodal mass is Σ_{two-phase p} W_pa m_p, × 2^k on a B-spline face layer;
  - `pressureCutoff` defaults to 5e-2;
  - m_max = max over all particles of m_p.
- **The weights.** Each two-phase particle's S^p and ∇S^p are renormalised over its p-active nodes
  by B9's rule, σ = Σ W, γ = Σ ∇W. Nothing is done when every kept node is p-active.
- **Every other node's p** is fixed at 0 by the host SP; such a node is absent from the projection
  and from G2P.
- **The T8 stray-p row**, a tail of 1e-8 giving p = −6.7e3 at step 2, becomes a pass.

### Route S as built (232l-C, 2026-11-18; implementation notes, nothing frozen moves)

- `Engine::advanceUP` (`LadrunoMPMRouteS.h`) runs every mode-F step of a model with a two-phase body,
  a plate, a `pfix` or a traction; any other model keeps `advance` bit for bit (T10).
- The supports are mode A's per-component folded ones (`upSupport`: P9, P10). M⁻¹ is §7a Q4's rule,
  the plate unknown included. The u side keeps mode F's cutoff and no renormalisation.
- p-activity, S^p and P5a are the rulings of mode A: the two-phase nodal mass against
  `pressureCutoff`, raw-weight renormalisation with folded-away entries retained (bb1421559), and
  the probe-guarded consistent projection. Particles carry p only (no ṗ, p̈ in mode F).
- The pressure solve: Jacobi-PCG on the free p-active rows, ||r|| ≤ 1e-13 ||b||, 4 n + 1000
  iterations, warm start p̂; drained p₀ on the right side. Non-convergence refuses the step.
- `mpm info gridPressure` returns [id, p, p̂] per p-active node of the last accepted step
  (T0 iii). `mpm criticalTimeStep drainedPencil | undrainedPencil`: 200 power iterations on
  M_L⁻¹K at the committed state; the undrained one adds α²Q̄ m mᵀ at two-phase particles.
- **Proposed amendment A7 (not applied silently):** route S refuses fix 2's affine image term
  (an impervious B-spline face with a normal fluid drive) in this build. No gated route-S case
  (T0 iii, T5, T6, T7) has one. Carrying π and the q + π invariant into mode F is offered as a
  follow-up.

## 8. Proposed amendments and objections (none applied)

- **A1** (U4): particles carry p̈ as well as p, ṗ (§1, P2).
- **A2** (T1): `mpm kinematics small` is new code that T1 presupposes; hypo stays the default.
- **A3** (U6): route S is (b)'s structure on u-p equations; Berkeley's own two-velocity,
  incompressible scheme is out of U1 and is cited, not ported.
- **A4** (U8, T3, T6): one mechanism, `mpm plate`, loads Terzaghi and Mandel; on bspline2 it rests
  on two new folds (P9's per-component map, P10's affine image). Each gets a patch check (a drained
  elastic block under the plate: uniform σ′_yy, u_y = U at every face node, rel 1e-10) before any
  T3/T6 number; if a fold fails it, that basis's Mandel legs become reported with the diagnosis.
- **A5** (T5): mode F has no Rayleigh; B4's 0.05M + 0.02K applies to route E only, a stated gap.
- **A6** (T7): the yardstick is the discrete drained pencil (`drainedPencil`), not h/c_d, whose
  2D bilinear value is not the stability limit; on the 1-cell GIMP column (all x fixed) they agree.
