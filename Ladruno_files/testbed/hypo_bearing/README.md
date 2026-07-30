# hypo_bearing — ADR-79 bearing-backbone campaign

Answers the ADR-78 §1 question with the now-complete geometry ladder
(`linear ⊂ corot ⊂ hypo ⊂ hypo -kozenyCarman`, ADR 79): does genuine large
strain (+ porosity/permeability evolution) bend the PDMY footing backbone
toward a bearing limit point, or reduce the reported ~9.3×-Vesic
over-hardening?

**Results and the cross-method verdict: `RESULTS.md`.**

## The model

A square **B = 2 m** smooth footing at the centre of the top face of a graded
all-hex box, saturated PDMY sand (medium-loose, φ = 33°), pushed
displacement-controlled to s/B = 0.15.

| | |
|---|---|
| domain | 20 × 20 m plan × 8 m deep = **4.5 B side clearance, 4 B depth** |
| mesh | 2816 `LadrunoUP` H8, 3468 nodes, 13 872 u-p DOF; 0.5 m hexes under and near the footing, graded outward (r ≈ 1.44 in plan, 1.35 in depth) |
| footing | the central 4 × 4-element patch (25 driven nodes), smooth (only u_z driven) |
| surcharge | 10 kPa over the **whole** top face, tributary-weighted, applied with gravity |
| staging | elastic gravity → `updateMaterialStage 1` → plastic settle → push |
| legs | ladder: `linear` (base rung), `corot`, `hypo`, `hypo -kozenyCarman`; locking: `corot_std`, `corot_bbar` (+ `_und` undrained pair) |

Normalization is Vesic **with** the surcharge term,
`q_ult = q0·Nq·sq + 0.5·γ'·B·Nγ·sγ` = 637.5 kPa (430.4 + 207.1).

## Scoping findings

Findings 1–4 were measured on the original uniform-box scaffold (2026-07-28);
findings 5–7 were measured while bringing the graded mesh up; finding 8 came
from the deep-push run (2026-07-29); finding 9 from scoping the undrained
locking pair (2026-07-30).

1. **Displacement control is mandatory** (from ADR-79 P3): a force-controlled
   push on stage-1 PDMY diverges from the first increment — near-surface GPs at
   ~zero confinement make the tangent nearly singular and the first Newton
   iterate unbounded.

   *Amended 2026-07-29 — the `-geom linear` half of this finding was too
   strong.* The original wording ("grinds without converging, excluded from the
   ladder") was measured with the pre-adaptive machinery: 2.5 mm increments
   driven by plain Newton, which we now know fails for **every** geometry
   method (finding 6). Re-run with the adaptive increment + KrylovNewton,
   `linear` converges perfectly well and is in fact the *best*-conditioned leg
   (0.15 retries/step vs 0.16 for hypo) — it is a usable base rung and the
   ladder is now measured against it. But the original verb was still the right
   one: it does eventually **grind**, just much later. The adaptive stepper
   moved the wall from push step 1 to s/B = 0.0469, where it dies by genuine
   DIVERGENCE (‖Δu‖ ~3e-2, residual 200–320) rather than the tolerance
   near-miss that ends `corot`. So: `linear` is admitted to the ladder as the
   base rung, with a hard ceiling at s/B ≈ 0.047.
2. **Surface-surcharge regularization is mandatory**, and it must cover the
   footing patch too. Without it the free lateral DOFs of zero-confinement
   surface nodes are near-singular under stage-1 PDMY. Re-measured on the
   graded mesh: an earlier build applied the surcharge only *beside* the
   footing (to keep q0 a purely adjacent overburden and avoid biasing the
   reaction) and every leg died on push step 1 — unbounded first Newton
   iterate, corot reporting inverted elements, then a KrylovNewton stall. The
   zero-confinement DOFs that need regularizing are the ones *under* the
   footing. The reaction bias this reintroduces is exactly
   `q0 × (footing-node tributary area)` and the runner subtracts it in closed
   form.
3. **Box confinement dominates every uniform-1 m-hex config that fits a quick
   run.** With roller sides at 0.5–1.5 B clearance the "footing" is an
   oedometer: q grows self-reinforcingly (PDMY G ∝ √p′ under rising
   confinement) to 27 MPa at s/B = 2.4 % — physically meaningless as bearing.
   This is what motivated the graded mesh (`build_mesh.py`); the real SFIM
   model is not in this repo.
4. ~~Per-step increments must stay ~2.5 mm~~ — **superseded by finding 6**.
   The 2.5 mm figure was calibrated on 1 m hexes and does not transfer.
5. **The engine must be asserted, not assumed.** An installed Ladruno wires a
   `.pth` into every venv that runs `import opensees` at interpreter startup,
   so a worktree `sys.path.insert` is a no-op and the campaign silently runs on
   the *installed* build — which predates ADR 78/79 and refuses
   `-geom corot|hypo`. Run with the base Python 3.12 (no Ladruno `.pth`); the
   runner carries a hard guard. See `LEDGER_quirks.md`.
6. **The converged push increment is a property of the mesh.** On this mesh
   2.5 mm does not converge for *any* geometry method — including
   `-geom linear`, which is what proves it is the model/mesh and not the
   geometry lane. It is **not** a tolerance artifact (Newton diverges at
   ‖Δu‖ ≈ 0.2 m; KrylovNewton stalls at ‖Δu‖ ≈ 3e-5 for `NormDispIncr`
   1e-8…1e-4, `RelativeNormDispIncr` and `EnergyIncr` alike), **not** a `dt`
   effect (2.5 / 25 / 250 s fail identically), and **not** the u-p
   stabilization or formulation (`-stab auto 0.10/0.25/0.50`, `-stab off`,
   `-formulation bbar` all fail identically at 0.25 mm).
7. **The step-size limit is a first-loading shock, and it is transient.** From
   the stage-1 gravity state 0.05 mm converges 12/12 while 0.10 mm fails; after
   a few small increments the same model accepts 0.4 mm (8×). So the push uses
   a 2-point linear ramp (settlement = rate × pseudo-time, increment carried by
   `dt`) with an **adaptive** increment — halve on failure, double after 5 good
   steps, truncate honestly at a floor. KrylovNewton is the *primary*
   algorithm, not a fallback; plain Newton diverged at every increment tested.
8. **How deep a rung can be pushed is itself an ordered result.** Each rung
   dies at a different penetration, in ladder order:
   `linear` 0.047 (divergence) < `corot` 0.060 (convergence floor) <
   `hypo` 0.150 (reached the target, still converging comfortably). Note the
   locking legs invert this within a rung: `corot_bbar` truncates at 0.036,
   *earlier* than `corot_std`'s 0.060, so relieving locking costs reach.
   Richer kinematics stay solvable longer on a
   mesh whose near-footing elements are being crushed — an argument for the
   rate-form UL lane that is independent of the backbone values themselves.
   Note `linear` does not taper into its wall: its last successful step was at
   the full 0.4 mm cap, and it then failed to complete the halving ladder in
   75 min because each diverged attempt costs ~10 min in PDMY's *serial*
   return mapping (~1 core busy vs ~5.3 for a healthy leg — a useful tell that
   a leg is diverging rather than merely slow).
9. **An undrained leg still needs a DRAINED initial state, and undrained-ness
   costs two knobs, not one.** Measured while scoping the `*_und` pair.
   *(a) Both knobs are required.* Undrained-ness is `T_v = c_v·t/B²`,
   `c_v = k·M/γ_w ≈ 2.28e4·k`. The drained legs sit at `T_v ≈ 8.5e4`, and
   permeability alone does not fix that: even clay-grade `k = 1e-8` only
   reaches `T_v = 8.5` — still drained — because the adaptive push spends
   ~1.5e5 s of pseudo-time getting to s/B = 0.06. *(b) The rate has a floor.*
   `V_s = √(G/ρ) = 166 m/s`, so a shear wave crosses B in 0.0121 s;
   `T_PUSH/1000` puts `dt` at the smallest increment at 0.005 s = 0.4× transit,
   and that leg measures inertial stiffening, not undrained locking.
   `T_PUSH/100` keeps `dt ∈ [0.05, 3.2] s` = 4.1–265× transit. So
   `tscale = 100` (all the rate the wave floor allows) plus `k = 1e-9` for the
   rest ⇒ `T_v = 0.0137`. *(c) Gravity must still consolidate.* At `k = 1e-9`
   the stock 500 s gravity `dt` leaves the GRAVITY stage undrained too: pore
   pressure carries the overburden, effective stress stays ~0, and
   `updateMaterialStage 1` then hands PDMY a near-zero-confinement tangent that
   will not settle (`plastic settle failed`) — finding 2's near-singular
   surface DOFs, reached through drainage instead of a missing surcharge. The
   fix is standard practice: consolidate drained, then shear undrained. Gravity
   `dt` is per-leg (`grav_dt`, 5e4 s undrained ⇒ `T_v = 5.7`), which brings the
   probe datum to p0 = 9.851 kPa against a hydrostatic 9.810. *(d) Verify the
   regime, don't assume it.* The runner records EXCESS pore pressure at an
   interior probe (footing centre, B/2 deep) over the consolidated datum —
   absolute p is useless here, since the probe sits 1 m down where hydrostatic
   alone is 9.81 kPa and raw `p/q` reads ~0.6 at ANY permeability. Measured
   `Δp_exc/Δq`: **0.000 drained** (0.003 kPa over 66.9) vs **0.360 undrained**
   (6.98 over 19.4).

## Files

- `build_mesh.py` — graded all-hex mesh builder (apeGmsh; 3 × 3 × 2 fragmented
  transfinite blocks, recombined). Writes `bearing_mesh.npz` (nodes, hexes,
  surcharge tributary areas, node sets). Run it with the **apeGmsh** env
  (`opensees_env`), which has gmsh + apeGmsh.
- `bearing_mesh.npz` — the cached mesh, i.e. the exact input to the CSVs here.
- `bearing_backbone.py` — the runner (staging idiom, surcharge, displacement
  control, adaptive increment, truncation-honest summary, incremental CSV).
  Run it with the **base Python 3.12** (finding 5).
- `backbone_<leg>.csv` — raw backbones (`linear` / `corot` / `hypo` /
  `hypo_kc`, plus the locking legs below), written incrementally. A run capped
  with `ADR79_MAXSTEP` is a smoke and writes `backbone_<leg>__smoke.csv`
  instead, and any leg refuses to open a CSV another process touched in the
  last 180 s (`ADR79_FORCE=1` overrides). Both guards exist because a smoke of
  a leg that was already running truncated the live file under it: the live
  process kept writing at its old offset, leaving the smoke's rows, a block of
  NUL padding, and only the tail of the real run.
- `leg_<leg>.log` — per-leg run logs.
- `RESULTS.md` — the write-up and the verdict.

### Leg groups

- `all` = the geometry ladder (`linear`, `corot`, `hypo`, `hypo_kc`).
- `pair` = `corot_std` + `corot_bbar`, the LOCKING probe. Every ladder leg runs
  `-formulation std`, so the ladder's corot→hypo deltas are the geometric
  content of a volumetrically locked mesh. `-geom hypo` refuses bbar
  (`OPS_LadrunoUP.cpp`, ADR 79 P2), but corot composes with it unchanged
  (ADR 78 §3.1), so corot is where the locking lever is measurable today.
  `corot_std` re-runs the `std` baseline on whatever build is in `dist/bin`, so
  the ratio is engine-matched rather than mixing formulation with engine — and
  it leaves the committed `backbone_corot.csv` intact.
- `pair_und` = the same two legs undrained (finding 9), where bbar's
  near-incompressibility lever is expected to be largest.

```bash
# mesh (apeGmsh env), then the four legs (clean interpreter)
C:/Users/nmb/venv/opensees_env/Scripts/python.exe Ladruno_files/testbed/hypo_bearing/build_mesh.py
```

```bash
C:/Users/nmb/AppData/Local/Programs/Python/Python312/python.exe Ladruno_files/testbed/hypo_bearing/bearing_backbone.py all
```

```bash
C:/Users/nmb/AppData/Local/Programs/Python/Python312/python.exe Ladruno_files/testbed/hypo_bearing/bearing_backbone.py pair
```

```bash
C:/Users/nmb/AppData/Local/Programs/Python/Python312/python.exe Ladruno_files/testbed/hypo_bearing/bearing_backbone.py pair_und
```

A group runs its legs in sequence. To run them concurrently (the campaign did;
~3–4 h each either way), launch one process per leg name instead.
