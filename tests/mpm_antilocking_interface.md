# WP-232j — anti-locking in rate form (232j-C, round 1: the freeze candidate)

Status: **proposed** by 232j-C, 2026-11-17; binding once the orchestrator says "frozen". Authority:
ADR 232 (§4 W1–W17, §5.2, §5.3 + 232f amendments, §5.4, §6, §8 + the 232j-kickoff amendment,
§9 phase 4) and wp-briefs.md "232j" (K1–K6, gates L1–L8). Anchors are quoted from
`origin/wp/232j-mpm-antilocking` @ `5955546ad`; line numbers drift. `tests/mpm_reference.py` and
`tests/mpm_implicit_reference.py` were not read. Decisions are numbered (J*n*) for review.
Conventions: `tests/mpm_interface.md` §1 (row-major F, tension positive, engineering shear).

## 1. K1 — the scheme: a cell-centre F-bar on the step increment

- **J1 — equations (rate form, K2).** For particle p, with the cell c(p) that holds its
  step-start centre (`Grid::cellOf`, `LadrunoMPMCore.h:218`, per axis, clamped to n−1), and d = ndm:
  ```
  ΔF_p  = I + Σ_{a∈supp p} Δu_a ⊗ ∇S_a(x_p^n)                (today's increment, any basis)
  ΔF_c  = I + Σ_{b∈corners c} Δu_b ⊗ ∇N_b(x_c)               (Q1 field of the cell's 2^d corners at its centre x_c)
  ΔF̄_p  = (det ΔF_c / det ΔF_p)^{1/d} ΔF_p ,   F_p^{n+1} = ΔF̄_p F_p^n   (the only F carried)
  ```
  Δu = Δt v* in mode F, the trial grid displacement in mode A. ∇N_b(x_c) = s_b/(2h) per axis in 2D
  and s_b/(4h) in 3D (s_b = ±1 for the corner's side) — the same constants for all four bases.
  Plane strain: 2×2 increments, d = 2, the embedded ΔF̄_zz = 1 by construction (L1c). The kernel
  entry is unchanged: `hypoIncrementFromF(F_n, ΔF̄ F_n, …)`; R̄₁ = polar(ΔF̄F_n) = polar(ΔF F_n)
  (a positive scalar does not change the polar factor). `det ≤ 0` on either side is a kernel
  refusal (today's −1 / −2 path). Total F is never averaged (K2).
- **J2 — averaging domain and weights.** The domain is the background cell; the "average" is the
  projection of the volumetric increment onto cell-wise constants, the constant being the
  cell-centre value of the Q1 interpolant — de Souza Neto et al. (1996)'s centroid F̄, applied to
  the background cell. No particle weights: every particle of a cell gets the same det ΔF_c. In
  2D the Q1 map's Jacobian is linear in (ξ, η), so det ΔF_c is exactly the cell's area ratio; in
  3D it is the one-point value. Coombs et al. (2018) introduced the incremental F-bar to MPM
  with a cell-constant dilatation (secondary sources: Xie et al. 2023, eq. 15 attributes the
  incremental form to them; Telikicherla & Moutsanidis 2022: "projection of the dilatational …
  part … onto a space of constant functions"); their exact averaging formula was **not verified**
  (paywalled; the Durham repository returned 403), so this file claims only that lineage.
- **J3 — per basis the corners are inside the support.** For x_p in cell [x_i, x_{i+1}] on an
  axis, both corners are at distance ≤ h: GIMP keeps nodes with |r| < h + l_p
  (`gimpFactor`, l_p > 0); `bspline2` / `bspline3` reach 1.5h / 2h with S(h) = 1/8 / 1/6 > 0;
  `linear` emits exactly the cell's two nodes, even at S = 0 (`if (b == BASIS_LINEAR) {`, `:250`,
  no drop). So **corners(c) ⊆ supp(p) before any cutoff, for all four bases**: J1 adds no node.
  For `bspline2` the spline field's own centre gradient equals the Q1 one exactly (dS(±h/2) = ∓1/h,
  S(±h/2) = 1/2, S(3h/2) = 0); for `bspline3` it differs at O(h²) — Q1 is used for every basis so
  the rule is one rule. Images (K4 of 232i) fold ghosts onto real nodes; corners are real nodes.
- **J4 — why L1(a) holds.** Every basis reproduces affine fields on the extended lattice
  (Σ∇S = 0, Σ x_a⊗∇S_a = I; 232i G1), and so does Q1. An affine Δu_a = A x_a + b gives
  ΔF_p = ΔF_c = I + A, so the factor is 1 to rounding: |ΔF̄ − ΔF| ~ 1e-16 against |ΔF − I| ~ |A|.
  L1 uses |A| ≈ 1e-2 (rel ≈ 1e-14 ≤ 1e-12). L1b: A = Q − I (30°) gives ΔF_c = ΔF_p = Q, so ΔF̄ = Q
  to rounding and the kernel's `rigid_rotation_*` objectivity applies at its frozen tolerance.
- **J5 — cells and particles.** *Changes cell:* the cell is that of x_p^n, like S and ∇S; a cell
  carries no history (rate form), so a crossing carries nothing. *One particle:* the same rule —
  the value comes from the grid, not from particle counts (an average over particles would leave
  a one-particle cell locked). *Empty cell:* nothing to compute. *A corner not active* (mode F
  `gact = 0`; mode A dropped by the cutoff / `renormalize`, excluded, or inactive in `ghost`
  mode): every particle of that cell **falls back to ΔF̄ = ΔF** (no averaging). This is a free
  surface's outer cells; it is the identity there, so L1 and L6 are untouched.
- **J6 — rejected candidates** (K4 is decisive):

  | Candidate | Coupling of particle p | 3D node set (gimp / bspline) | Fate |
  |---|---|---|---|
  | J1 cell-centre F-bar | own support | unchanged: 8–27 / 27 / 64 | **chosen** |
  | volume-weighted mean of the cell's particles' ΔJ | union of the cell's supports | ≤ 4/axis → 64 nodes, 192 DOF | rejected: 192² × 8 B = 295 kB own tangent per particle-element; one-particle cell locked |
  | grid projection ΔJ_I = Σ N m ΔJ / m_I, back by N (Zhao et al. 2022; Xie et al. 2023) | 2-ring | ≤ 7/axis → 343 nodes, 1029 DOF | rejected: unboundable; its users skip the linearisation (Xie et al.: Broyden) |
  | reduced-degree B-bar/F-bar projection (Elguedj et al. 2008; Telikicherla & Moutsanidis 2022) | projection support | larger still, and splines only | rejected: not all four bases |

## 2. Mode F — where and what it costs

In `Engine::advance`, step 5, after the second projection and its conditions
(`for (size_t s = 0; s < slave.size(); s++) reflectVec(gvs, s);`, `:855`) and inside the particle
loop that forms `L` (`:861`), between `dF` and the kernel call (`:882`) — i.e. after the grid update,
before `mats.setTrialStrain` (`:893`). ΔF_c reads `gvs` and `gact` at the 2^d corners. Volume
follows: `tV[p] = J1 * P.V0[p];` (`:907`) now holds det(ΔF̄F_n)V0, and the next P2G
(`gfi[nd * i + a] -= Vp * s;`, `:744`) uses that carried volume and σ̄. Cost per particle-step:
2^d·d² multiply-adds, two determinants, one `sqrt`/`cbrt`, d² scalings (≈ 100 flops). Estimate
2–8 % of an elastic step; measured in R (wall per step, on vs off). Same thread-determinism as
today (per-particle, no reduction).

## 3. Mode A (K4) — no new coupling, no new element

- **J7 — dependence.** f_p depends on Δu over supp(p) ∪ corners(c) = supp(p) (J3). The element's
  node set, `numDOF`, re-pointing and Domain churn are **unchanged**; no per-cell element.

  | basis | 2D nodes / DOF | 3D nodes / DOF | 3D vs `MAX_NUM_DOF` = 64 (`FE_Element.cpp:48`) |
  |---|---|---|---|
  | linear | 4 / 8 | 8 / 24 | pooled |
  | gimp | 4–9 / 8–18 | 8–27 / 24–81 | 27 nodes own a tangent (`:188` vs `:210`, #1098) — as today |
  | bspline2 | 9 / 18 | 27 / 81 | own tangent — as today (#1098) |
  | bspline3 | 16 / 32 | 64 / 192 | own tangent; element `MAX_DOF = 192` (`LadrunoMPMPoint.cpp:54`) |

- **J8 — bytes (M0b M3 method).** New fixed members, present on and off: `int corner[8]` (slot in
  `connect`, −1 = fallback), `double Fci[9]` (ΔF_c⁻¹), `double Jc, hc` ≈ 136 B, +0.8 % of the
  17.8 kB measured (L6: ≤ 5 % off). No per-element matrix: K_F is added in the static scratch.
- **J9 — residual.** Today's lines with ΔF̄: `Frel` (`:240`) stays the actual ΔF and still gives
  g_a = ∇S_a·ΔF⁻¹ (`inv3(Frel, Fi)`, `:282`) and the A7 trial centre; the kernel takes ΔF̄F_n;
  `vol = J * PP.V0[p];` (`:278`) becomes w̄ = det(ΔF̄F_n)V0, the **carried** volume. Hence
  f_a = w̄ σ̄ g_a. Choosing the carried volume (not de Souza Neto's actual J) keeps today's property
  that a converged step's internal force equals the next step's at Δu = 0, and is mode F's
  quadrature weight: one scheme in both modes. G2P (`:425-426`) then stores F̄ and w̄ unchanged.
- **J10 — consistent tangent, geometric term included.** With ℓ = Σ δu_b⊗g_b, the hypo model
  that makes today's K exact, δ(J̄σ̄) = J̄(c:d̄ + ℓ̄σ̄ + σ̄ℓ̄ᵀ) with ℓ̄ = ℓ + βI_d,
  β = (1/d)(tr ℓ_c − tr ℓ), gives
  ```
  K = w̄ [ B_aᵀ c̄ B_b + (g_a·σ̄·g_b) I ]  (today's assembly at σ̄, c̄ = pushTangent(R̄₁, D(F̄)))
    + h dᵀ ,   h_ai = (w̄/d) Σ_j [ c̄:I_d + 2σ̄ ]_ij g_aj ,   d_bk = ĝ_bk − g_bk ,   ĝ_b = ∇N_b(x_c)·ΔF_c⁻¹ (corners), 0 elsewhere
  ```
  (c̄:I_d)_r = Σ_{q<d} c̄[r][q] in Voigt. K_F is **rank one** (O(n²), n ≤ 192) and unsymmetric even
  for elastic materials (as in de Souza Neto et al. 1996). B8's per-component path uses g′_{a,i}
  and d_bk = ĝ_bk − g′_{b,k,k}. `getInitialStiff` adds the same term at Δu = 0, σ̄ = 0, initial D.
  Fallback particles add nothing. #1096's B-build is not rewritten (K_F is a separate loop).
- **J11 — re-pointing and empty nodes.** `prepareStep` computes c(p) at x_p^n, maps its 2^d corner
  ext ids to slots in the kept tag list; a corner absent from it (renormalize / exclude) or
  inactive (ghost) sets slot −1 for the whole particle (J5). `repoint(…)` gains `const int
  *corner` (null when off: the shipped path, bit for bit).
- **J12 — cut and retry.** Nothing new to undo: ΔF_c is recomputed from the trial Δu in every
  `update()`; G2P (batched, D1.9) writes F̄ only at commit. A retry re-prepares from the
  committed particles. #1095 is not touched (and still a hazard for the ladder, §5).
- **J13 — restart pack.** `MPM_PACK_VERSION` 2 → 3 (`:1145`), the setting appended after the body
  force (`k.putv(eng.bodyForce);`, `:1184`). `unpackState` accepts 3, and 2 as **off**; 1 stays
  refused as today (#1103 item 2 is not in scope and not regressed). The element pack is
  unchanged (corners are rebuilt at the next re-point; `hdr(8)`, `:769`, stays 0).

## 4. K3 — syntax, default, refusals

- `mpm antiLocking off | fbar` — default `off`. Fixed for an epoch like `mpm mode`
  ("the mode is fixed at the first step of an epoch; `mpm reset` first", `LadrunoMPMImplicit.cpp:174`).
  Kept by `mpm reset`; gone with `wipe`. Query: `mpm info antiLocking` → 0 / 1; `Print` names it.
- **Off is bit for bit:** every new branch is guarded by the setting (mode F) or by a null corner
  list (mode A); no arithmetic on the off path changes.
- **Refused, loudly:** an unknown keyword; a change after the first step of an epoch; `fbar` in
  mode A with a symmetric-storage SOE — `BandSPD` (3), `ProfileSPD` (4), `SymSparse` (11),
  `ProfileSPDGather` (14) (`classTags.h:1188-1199`) and PARDISO in its symmetric mode (flag quoted
  in round 2) — because K_F is unsymmetric and those would give a wrong answer silently. The SOE
  is reachable: `theSOE` (`commands.cpp:653`), `cmds->getSOE()` (`OpenSeesCommands.h:87`); the
  bridge gains it. No basis, mode or face is refused: J1–J5 cover all.
- `mpm debug tangentCheck <scale> <eps> [-antiLocking off|fbar]` evaluates either setting at the
  same committed state (L2); `mpm debug affineCheck <A: d²> <b: d>` (both modes, the
  `tangentCheck` arming pattern in mode A) puts Δu_a = A x_a + b on every support node, returns
  [max_p |ΔF̄−ΔF|/|ΔF−I|, max_p |Δε̂|, #averaged, #fallback] (L1).

## 5. K6 — the footing ladder on mode A

The harness owns it (232j-V's `ladder` key), per accepted step: `mpm cut 2 0` once (no internal
cuts); for each try k: `algorithm <alg>`, `test NormUnbalance <tol> <it>`, `mpm analyze 1 h`.
A failed `mpm analyze` with 0 cuts reverts to the last accepted step and aborts nothing
(`if (++cuts > maxCutCount) {`, `:1050`, rc −1 without `aborted`), so the next try re-prepares
cleanly. Try ≥ 3 accepted = relaxed. All tries fail → h/2, regrown ×2 per accepted sub-step to dt;
the quasi-static driver sets Δλ = h (`ladrunoSetIncrement(h, h, h, -1.0)`, `:1031`) and gridvel
gives Δu = v·h. Counted per run: failed attempts, relaxed steps, halvings. Same in both engines
(a Tcl proc). System: UmfPack (already the harness's, unsymmetric).

## 6. K5 — contracts unchanged

The setting changes no field's meaning, shape or units, and is fixed per epoch (§4), so a
stage is never mixed. C2-GRID records the basis because a reader needs it for grid metadata;
nothing here is needed to read results. Provenance lives in the deck and the results JSON. An
optional additive attribute is a T follow-up only if the owner asks.

## 7. Harness plan

`tests/mpm_antilocking_cases.txt` (every case is **on** unless it says `antilocking off`, so
232j-V's L3 block needs no edit): **L1** `affine_{F,A}_{basis}_{2d,3d}` (3^d cells, ppc 2, body =
real grid; `masscut 1e-12` so no support is truncated), deviatoric, volumetric and 30° cases;
**L2** `column_selfweight_small_qs` undeformed and at 0.2 % strain, on vs off; **L4**
`footing_mildcone_{r1,r2,r3}_on` (bands [127.29, 155.58], [117.60, 143.73], [110.28, 134.78];
`reach_t 110`; q/exact, plateau, s at q_max, ladder counts reported); **L5** G2 / G4 on (rows of
`mpm_bspline_cases.txt`); **L7** `switch_j2_r1_{on,off}` (§8 item 3).

| test file (zone_a) | rows | Zone-A cost (est.) |
|---|---|---|
| `test_ladruno_mpm_antilocking_kinematics.py` | L1 g++ (`tests/mpm_antilocking_core_check.cpp`) + engines; L8 (info, refusals, restart v2→v3, Tcl = Python on one L1 case); L6 off spot rows | ~30 s |
| `test_ladruno_mpm_antilocking_tangent.py` | L2 (ε 1e-8, 1e-9, 1e-10) | ~5 s |
| `test_ladruno_mpm_antilocking_switch.py` | L7 switch oracle, on and off | 2 × footing r1, 60 steps (232f runs `footing_mildcone_r1`, 110 steps, in Zone-A today); measured round 2 |
| `test_ladruno_mpm_antilocking_gates.py` | L3, L4, L5 full ladders | skipped unless `LADRUNO_MPM_GATES=1` (warrant run) |
| `test_ladruno_mpm_antilocking_results.py` | asserts the gates on `tests/data/mpm/antilocking/*.json` (ladrunoBuild, ladder counts) | < 1 s |

`test_footing_gate` (`test_ladruno_mpm_implicit_gates.py`) becomes the L4 test; the old band's
verdict is printed, never asserted. L5's Zone-A smoke: rungs 8–32, non-increasing only (≈ 7 s).

## 8. Proposed amendments (none applied)

1. **L1(a) scope:** "every particle whose support is not truncated by the empty-node treatment".
   A truncated (renormalised) support does not reproduce affine fields, so its own ΔF is not
   I + A while ΔF_c is; the cases avoid truncation and the truncated result is reported.
2. **L8:** `mpm info summary` carries no settings (`[mass, momentum, max |v|, KE]`); the round trip
   is `mpm info antiLocking` (the 232i G7 precedent).
3. **L7 switch oracle band:** L3's closed-form band is an r3 gate and even FE B-bar misses it at r1
   (570.10 = 1.109×), so at r1 the oracle uses L3's fallback form: on within ±10 % of 570.10 →
   [513.09, 627.11]; off outside it (FE std 650.72 is). Fixed here, before any MPM number.
4. **Citations (K1):** "Moutsanidis et al. 2020" is IGA-MPM (CMAME 372, 113346) — verified, but it
   is not the reduced-degree locking projection; that is Elguedj et al. 2008 and, in MPM,
   Telikicherla & Moutsanidis 2022.
5. **New refusal (K3):** symmetric-storage SOEs with `fbar` in mode A.

## 9. Sources (literature only; each verified 2026-11-17)

de Souza Neto, Perić, Dutko & Owen (1996), IJSS 33(20–22):3277–3296, doi 10.1016/0020-7683(95)00259-6.
Coombs, Charlton, Cortis & Augarde (2018), CMAME 333:1–21, doi 10.1016/j.cma.2018.01.010.
Xie, Navas & López-Querol (2023), IJNAMG 47:2741–2761, doi 10.1002/nag.3599.
Zhao, Jiang & Choo (2022), arXiv:2209.02466. Elguedj, Bazilevs, Calo & Hughes (2008), CMAME
197:2732–2762. Telikicherla & Moutsanidis (2022), CMAME 395:114985. Hughes (1980), IJNME
15(9):1413–1418. Moutsanidis, Long & Bazilevs (2020), CMAME 372:113346.
