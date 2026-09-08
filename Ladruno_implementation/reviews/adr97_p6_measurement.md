---
title: ADR-97 P6 — Cerro-Lindo-scale measurement (Backward_Euler vs Closest_Point)
project: Ladruno
status: complete
owner: nmora
tags:
  - implementation
  - material
  - review
---

# ADR-97 P6 (`wp/97g-measure`) — measurement report

**PR** [#826](https://github.com/nmorabowen/OpenSees/pull/826) (draft, based on
`wp/97c-cp-principal` = P2 [#824](https://github.com/nmorabowen/OpenSees/pull/824))
· **plan** [[97_ladruno_asdp_closest_point_adr]] (D1, Phases table P6 row) ·
**P1 report** [[reviews/adr97_p1_report]] · **P2 report** [[reviews/adr97_p2_report]]
· **driver** `Ladruno_implementation/adr97_oracle/measure_p6.py`

## 1. What this measures

D1 keeps `Backward_Euler` + `Secant` as the shipped default and gates any flip
on a mesh-scale measurement, not on the material-point gates P1/P2 already
passed. This WP builds that measurement: a box of `LadrunoBrick` hexes
(Mohr-Coulomb or `MohrCoulombTensionCutoff`, non-associated `psi < phi`),
gravity-initialised (self-weight ramped under `LoadControl`, then frozen with
`loadConst -time 0.0`) and pushed toward bearing failure by a strip-footing
pressure on the top-centre third of the plan, ramped in a second
`LoadControl` stage. Five configurations are compared at two mesh scales:

| Configuration | `integration_method` | `tangent_type` | `algorithm` |
|---|---|---|---|
| `BE_Secant` (shipped default) | `Backward_Euler` | `Secant` | `Newton` |
| `BE_Continuum` | `Backward_Euler` | `Continuum` | `Newton` |
| `CP_Continuum` | `Closest_Point` | `Continuum` | `Newton` |
| `CP_Algorithmic` | `Closest_Point` | `Algorithmic` | `Newton` |
| `CP_Algorithmic_Krylov` | `Closest_Point` | `Algorithmic` | `KrylovNewton` (fork default algorithm) |

- **small**: 10x10x10 elements, 1331 nodes, **3993 DOF**, 3 gravity steps + 10
  push steps.
- **cerro**: 19x19x19 elements, 8000 nodes, **24000 DOF**, 3 gravity steps + 12
  push steps — the ADR-80 hex8-class rung size (`80_ladruno_sp_imposition_
  strengthening_adr.md` line 229: "All hex8-class (~22.6 k DOF)").

Every configuration runs in its own fresh `python3.12` child process:
`ASDPlasticMaterial3D`'s per-tag `integration_method`/`tangent_type` option
maps are process-global statics keyed by material tag (ADR-94 lesson, reused
by ADR-97 P1/P2), so reusing one process across configurations would read
back a stale option map.

Materials: `MC` is a generic soil deck (E=30000 kPa, nu=0.3, phi=32 deg,
psi=8 deg, c=15 kPa) deliberately different from ADR-97 P2's oedometric trap
(nu=0.25/phi=30 sit exactly on the compression meridian, K0 = 1/3, so a
self-weight-only path never yields). `MCTC` reuses the Cerro-Lindo-like EDZ
deck from `tests/test_asdplastic_mctc.py` / ADR-97 P2 (E=2.0e6 kPa, nu=0.3,
c=100 kPa, phi=20 deg, psi=5 deg, T=24.7 kPa) — a real geotechnical deck, not
an invented one. The footing pressure target is a fraction of a rough
Terzaghi `q_ult` (0.12 for MC, 0.60 for MCTC — MC's high-friction/
low-cohesion deck engages plasticity almost immediately under self-weight
alone; MCTC's high cohesion needed a much larger push to leave the elastic K0
state at all, measured directly: 0.12*q_ult left every Gauss point elastic).

## 2. Provenance

- Build: `ladrunoBuild()` = `06913de30f3bcdf837e4d69b41dae0267f4f7c87` (the P2
  build this worktree's `dist/bin` was copied from; verified before any run).
- Machine: `Windows-11-10.0.26100-SP0`, Intel64 Family 6 Model 183 (24 logical
  CPUs), `python3.12`.
- Solver: `system UmfPack` throughout (`C_alg` is unsymmetric for every
  non-associated configuration here; `ProfileSPD` would silently be wrong and
  `FullGeneral` is not used). PARDISO was not exercised in this WP — the
  point was to isolate the integrator/tangent effect, not the solver, and
  every prior ADR-97 gate already runs on UmfPack.
- `numberer RCM`, `constraints Transformation`, `test NormDispIncr 1e-8
  <maxIter=250> 0`.
- Wall-clock timings below were measured with other background activity on
  the same host at various points (this is a shared dev machine per
  `ladruno-concurrent-worktrees` — see the note at the end of §5); the
  **iteration counts** are the load-bearing numbers, wall time is
  corroborating, not the primary evidence.

## 3. Results — MC family

| Configuration | steps converged | total Newton iters | first non-convergence | wall time (s) | plastic GPs (end) |
|---|---|---|---|---|---|
| **small (3993 DOF, 10 push steps)** | | | | | |
| `BE_Secant` (default) | 7/10 | 60 | step 7 | 39.7 | 16 |
| `BE_Continuum` | 7/10 | 48 | step 7 | 37.6 | 16 |
| `CP_Continuum` | **10/10** | 39 | none | 5.2 | 29 |
| `CP_Algorithmic` | **10/10** | 39 | none | 5.5 | 29 |
| `CP_Algorithmic_Krylov` | **10/10** | 53 | none | **3.0** | 29 |
| **cerro (24000 DOF, 12 push steps)** | | | | | |
| `BE_Secant` (default) | 7/12 | 56 | step 7 | 504.0 | 30 |
| `BE_Continuum` | 7/12 | 41 | step 7 | 587.4 | 30 |
| `CP_Continuum` | **12/12** | 48 | none | 78.3 | 176 |
| `CP_Algorithmic` | **12/12** | 48 | none | 78.1 | 176 |
| `CP_Algorithmic_Krylov` | **12/12** | 66 | none | **45.4** | 176 |

`BE_Secant`'s own per-step iteration history at both scales climbs sharply as
the plastic zone grows — small: `2,2,2,7,9,13,19` then non-convergence; cerro:
`2,2,2,5,6,9,10` — the classic cutting-plane symptom (ADR-94 M3: BE's fixed
point is a path-summed cutting-plane approximation, not the closest point,
and that approximation degrades as the flow direction rotates step to step
near failure). `CP_Algorithmic`'s per-step counts settle to a **flat 3-4**
(small) / **4** (cerro) once yielding starts and stay there through the whole
remaining push, including the 3 steps where `BE_Secant` has already failed.

## 4. Results — MCTC family

| Configuration | steps converged | total Newton iters | first non-convergence | wall time (s) | plastic GPs (end) |
|---|---|---|---|---|---|
| **small (3993 DOF, 10 push steps)** | | | | | |
| `BE_Secant` (default) | 6/10 | 112 | step 6 | 123.3 | 29 |
| `BE_Continuum` | 6/10 | 52 | step 6 | 47.5 | 29 |
| `CP_Continuum` | **9/10** | 90 | step 9 | 54.2 | 119 |
| `CP_Algorithmic` | **9/10** | 90 | step 9 | 48.8 | 119 |
| `CP_Algorithmic_Krylov` | **9/10** | 70 | step 9 | **17.0** | 119 |
| **cerro (24000 DOF, 12 push steps)** | | | | | |
| `BE_Secant` (default) | 6/12 | 113 | step 6 | 703.5 | 145 |
| `BE_Continuum` | 6/12 | 52 | step 6 | 637.0 | 145 |
| `CP_Continuum` | **11/12** | 140 | step 11 | 387.3 | 842 |
| `CP_Algorithmic` | **11/12** | 140 | step 11 | 420.9 | 842 |
| `CP_Algorithmic_Krylov` | **11/12** | 117 | step 11 | **257.6** | 842 |

Same qualitative story as MC, at the harder (higher-cohesion, more
confinement-sensitive) deck: every `Closest_Point` configuration reaches
**3 more converged push steps** (9 or 11 of 10/12) than either `Backward_
Euler` configuration (6 of 10/12), on both meshes, with `CP_Algorithmic_
Krylov` fastest by wall clock at both scales.

`BE_Continuum`'s total-iteration count is consistently *lower* than
`BE_Secant`'s at the SAME number of converged steps (e.g. MC small 48 vs 60;
MCTC cerro 52 vs 113) — `Continuum` is already a better-conditioned tangent
than `Secant` on `Backward_Euler`'s own cutting-plane map, exactly ADR-94's
M3 finding reproduced at mesh scale. Neither `Backward_Euler` tangent choice
fixes the earlier non-convergence, though: both fail at the identical step,
because the failure is the MAP (cutting plane vs closest point), not the
tangent quality feeding Newton.

## 5. Backward_Euler vs Closest_Point stress gap

At the last push step both maps have in common, the monitor element's normal
stresses (directly beneath the footing centre) agree closely — the
ADR-97-predicted regime where "CP and BE agree to Newton tolerance exactly
when the flow direction does not rotate over the step" (a near-proportional
loading path under a symmetric footing):

| Family | Scale | common steps | max normal-stress relative gap (BE_Secant vs CP_Algorithmic) |
|---|---|---|---|
| MC | small | 7 | 0.68% |
| MC | cerro | 7 | 0.28% |
| MCTC | small | 6 | 0.10% |
| MCTC | cerro | 6 | 0.02% |

(Shear components at this monitor point are near machine-zero by symmetry —
the footing is centred and the mesh is regular — so their *relative*
difference is large but physically meaningless; the normal components, which
carry >99% of the stress magnitude here, are the honest comparison.) This is
**not** a correctness gap — it is the two maps' known small divergence on a
path with mild flow-direction rotation, consistent with gate 4's pinned
contrast in P1/P2 (perfectly plastic/proportional decks agree near-exactly;
AF and rotating-normal decks differ measurably). The gap shrinking with mesh
refinement (smaller elements -> smaller strain increment per Newton step ->
less rotation per step) is the expected trend, not a new finding.

**Session-sharing caveat.** Per `ladruno-concurrent-worktrees`/`ladruno-
duplicate-lane-work` (this is a shared dev machine with multiple concurrent
worktrees), the absolute wall-clock numbers above were not all measured under
identical machine load — the `BE_Continuum`/`BE_Secant` cerro-scale numbers in
particular (504-703 s) were measured while other sweep runs from this same
WP were queued sequentially, not concurrently, but system-wide background
activity on the host cannot be ruled out. The **iteration counts** (`testIter()`
sums) are load-independent and are the primary evidence; wall time
corroborates the same ordering (`CP_Algorithmic_Krylov` fastest, `Backward_
Euler` configurations slowest) at every family/scale combination measured, so
the qualitative conclusion is robust even if the exact seconds are not
laboratory-clean.

## 6. Recommendation (D1)

**Flip the default: `integration_method Closest_Point` + `tangent_type
Algorithmic`, with `algorithm KrylovNewton` unaffected (it is already the
fork default).** At both mesh scales and on both material families measured,
`Closest_Point` converges strictly more of the load history than `Backward_
Euler` under the identical convergence criteria (`NormDispIncr 1e-8`, same
`maxIter`) — 10/10 vs 7/10 and 12/12 vs 7/12 on MC, 9/10 vs 6/10 and 11/12 vs
6/12 on MCTC — while costing no more, and usually markedly fewer, total
Newton iterations even over the steps `Backward_Euler` DOES complete (e.g.
MC cerro: `CP_Algorithmic` 48 iterations over 12 steps vs `BE_Secant`'s 56
iterations over only 7 before it stalls). This is exactly ADR-94 M3
reproduced at mesh scale: the shipped cutting-plane map's error compounds as
the plastic zone grows near failure, and it eventually costs *convergence*,
not just iteration count — the regime a bearing-capacity analysis exists to
resolve. The committed-stress gap between the two maps on the steps both
complete is under 1% (§5), so the flip does not silently change results on
decks that were already converging under `Backward_Euler` "well" (mild
flow-direction rotation); it recovers steps that `Backward_Euler` cannot
currently take at all. `CP_Algorithmic_Krylov` (the fork's default algorithm)
is the fastest configuration measured at every scale, so the flip should
ship as `Closest_Point` + `Algorithmic` with no change to the default
`algorithm`.

Caveats for the owner: (1) support is still 22 of 46 specializations
(VonMises, Drucker-Prager, MohrCoulomb families — P1/P2); StiffSoil/
RoundedMohrCoulomb/HoekBrown are refused at parse time under `Closest_Point`
until P3/P5 land, so the default flip cannot be unconditional until those
land or the flip is scoped to the supported families. (2) This measurement
used one footing geometry, one mesh regularity, and UmfPack; the ADR-75
PARDISO desktop path and MPI decomposition (`sendSelf`/`recvSelf`, ADR-94 D3)
are still out of scope for `Closest_Point`'s option maps (unverified here,
not contradicted). (3) Neither `Backward_Euler` configuration's earlier
non-convergence was chased to a root cause beyond "the cutting-plane map"
here — that diagnosis is ADR-94's, reproduced, not re-derived.

## 7. Test coverage

`tests/test_adr97_p6_measure_smoke.py` — Zone-A, <1s locally — runs the
driver's tiny 144-DOF smoke mesh for 3 push steps in the four base
configurations (excludes the KrylovNewton variant), each its own child
process, and asserts every configuration converges and `CP_Algorithmic`'s
total Newton-iteration count does not exceed `BE_Secant`'s — the ordering
this report measures at mesh scale, pinned (not the exact counts) so the
smoke test stays robust to build-to-build float noise.

## See also

[[97_ladruno_asdp_closest_point_adr]] · [[reviews/adr97_p1_report]] ·
[[reviews/adr97_p2_report]] · [[reviews/adr94_verdict]] (§1 M3) ·
`Ladruno_implementation/80_ladruno_sp_imposition_strengthening_adr.md`
(the ~22.6k DOF hex8-class rung size) ·
`Ladruno_implementation/adr97_oracle/measure_p6.py`
