# ADR-79 bearing campaign — session handoff (updated 2026-07-30, collapse study)

Status: **the campaign is complete and merged** (PRs #684, #686, #689, plus
#687/#688 on the concurrent locking line). The **collapse study — priority 1 of
the previous handoff — is done**; this note now carries its answer and the
priorities that survive it.

## Where everything lives

| what | where |
|---|---|
| runner, probes, figure scripts, CSVs | `Ladruno_files/testbed/hypo_bearing/` |
| geometry-ladder write-up | that folder's `RESULTS.md` |
| **collapse-study write-up** | that folder's **`COLLAPSE.md`** |
| fork decision record | ADR 79 **§8** (ladder) and **§9** (collapse) |
| traps | `LEDGER_quirks.md` (6 entries from this line) |
| **project-facing analysis** | TIMs vault notes **16** and **17**, `C:\Users\nmb\Dropbox\obsidian\Ladruño\TIMs\Reference Model\` — **both predate the correction below** |

Worktree `C:\Users\nmb\Documents\Github\OpenSees-hypo`, branch
`claude/adr79-bearing-deeper`, built `dist/bin` current.

## The answer, in four lines

1. **The geometry ladder is exonerated.** −5.4 % (corot) to +7.5 % (hypo)
   against a `-geom linear` base, no limit point out to s/B = 15 %, tangent flat
   at 6847 kPa/m. ADR-61 **OQ8 answered negatively**.
2. **The capacity anchor IS wrong — by ≈ 2.4×, not 52×.** The previous handoff
   said 52×, from Chen–Han plane-strain matching of the measured cone. That
   matching assumes **associated flow**, which this non-dilatant PDMY set
   (`dil1 = dil2 = 0`) violates. With Davis's ψ = 0 reduction the plane-strain
   equivalent is 38.87°, not 53.72°, and the square-footing anchor is
   **1525 kPa** — against the 637.5 kPa quoted, and against 26 711 kPa from the
   associated route.
3. **The measured collapse load agrees with the corrected anchor.** An
   elastic–perfectly-plastic Drucker–Prager on the *same measured cone*, same
   mesh, same footing gives **1140 kPa at the s/B = 10 % criterion** (1362 kPa
   extrapolated) = 0.75–0.92 of the Davis anchor and **1/23 of the Chen–Han
   prediction**. The machinery is validated to **1.0020 against the exact
   Prandtl–Reissner `q₀·N_q`**.
4. **So the over-strength is re-scaled, NOT dissolved.** The benchmark's
   3384 kPa at s/B = 15 % is **2.2× the correct anchor and ~3.0× the measured
   collapse load of its own cone**. The previous handoff's "the model sits at
   ~31 % of its own collapse load" is **withdrawn** — it is above it. Locking,
   large-strain embedment, boundary confinement and mesh coarseness still have
   to account for ~2.2×.

## What to do next, in priority order

1. **Re-anchor the vault ratios to 1525 kPa** (notes 16/17 quote the 52× route
   and the "31 % of capacity" reading; both are superseded). The recorded
   ratios themselves are fine as data — it is the interpretation that moves.
2. **A modelling decision for the project, sharpened.** A Lode-independent cone
   calibrated in triaxial compression is still the constitutive assumption in
   force, and it still over-predicts plane-strain strength — but by ~2.4× on a
   correctly-reduced anchor, not ~50×. Either calibrate φ for the regime of
   interest or use a Lode-dependent surface. Note the Chen–Han matching this
   model needs is **numerically unstable at this cone**: it has a pole at
   α = 1/√12 and the cone sits at 84.4 % of it, so ±2 % in α moves N_q by 55 %.
   That instability is an argument against the cone itself, not just against
   the formula.
3. **A refined mesh — now the binding constraint, and it needs REGULARIZATION,
   not more elements.** `COLLAPSE.md` tried refinement and it failed in a
   specific, diagnostic way: at 4 / 8 / 16 elements across B the
   non-associated Prandtl leg terminates at s/B = 0.0124 / 0.0030 / 0.0014 with
   the tangent at termination *rising* 26 % → 35 % → 63 %. A perfectly plastic
   non-associated frictional solid is past the Rudnicki–Rice localization
   threshold, so the continuum problem has lost ellipticity and band width is
   set by the element size. The next attempt should add a rate-dependent /
   viscoplastic regularization (or run it explicit-dynamic), not just refine.
   Everything quantitative in §9 rests on this gap. `-geom hypo` still refuses
   `-formulation bbar` (ADR 79 P2 reserved); `corot` composes with `bbar` today.
4. **The domain question at collapse — now known to be live.** At the full
   s/B = 0.15 the fully-mobilised zone **reaches the roller sides** (16 of the
   352 elements in the outermost column at m > 0.99; the base stays clear).
   `build_mesh_big.py` (14.5 B clearance, 10 B depth, 8064 hexes) agrees with
   the campaign mesh to 3 % out to s = 16 mm but has not reached the plateau —
   it is ~3× the cost per step and needs a long uncontended run. The sign is
   known (confinement inflates capacity), so finishing it can only widen the
   gap to PDMY's 3384 kPa.
5. **The interface**, once the contact `ndf == 3` guard is relaxed. Still the
   last untested kinematic candidate, and this benchmark cannot speak to it
   (smooth driven node patch, so an interface is absent by construction).

## Traps this line paid for — do not rediscover

All six are in `LEDGER_quirks.md`. The ones that cost the most time:

- **An installed Ladruno hijacks `import opensees`** via a site `.pth` that
  imports at interpreter startup, so a worktree `sys.path.insert` is a no-op.
  Use base Python 3.12 and **assert `opensees.__file__`**. A location guard is
  necessary but not sufficient — a worktree's `dist/bin` can be months behind
  its branch, so assert CAPABILITY too.
- **A frictional verification model can be 95 % mobilised by its own gravity
  state**, because the ELASTIC `K0 = ν/(1−ν)` sets the initial stress ratio.
  Compute `m = (1−K0)/(√3·α·(1+2K0))` and print it; `m > 0.8` is void. Raise ν
  for verification legs — a collapse load does not depend on elastic constants.
- **`DruckerPrager` cannot have pressure-dependent moduli AND plasticity**
  (`mElastFlag == 1` is elastic-only), and `updateMaterialStage` does not reach
  it (it answers to `materialState`). `ρ = √2·α`; perfect plasticity needs
  `H = 0`.
- **Consolidate-then-deviator with a `Linear` series applies the deviator at
  FULL amplitude on the first increment** — both patterns share the load
  factor. `loadConst("-time", 0.0)` between stages.
- **`NormDispIncr` is not evidence of equilibrium** (it reported convergence at
  σ_zz = +0.26 kPa tension under 100 kPa of compression). Use `NormUnbalance`.
- **A converged push increment is a property of the mesh**, not the problem.
  Drive with a ramp + adaptive increment; KrylovNewton is primary; and for
  perfectly plastic collapse add an ALGORITHM LADDER per increment
  (KrylovNewton → NewtonLineSearch → relaxed-tolerance KrylovNewton) before
  shrinking the step, flagging relaxed steps in the CSV.
- `remove('sp', node, dof)` does **not** match a constraint made by `fix`.
  UmfPack `setSize` fails at 4 DOFs (use `FullGeneral`). Field plots must be
  referenced to the gravity frame.
- **Two long parallel runs writing incremental CSVs**: re-running one leg's
  name truncates the live file under it. Smoke runs get their own filename, and
  a runner refuses to open an output modified in the last 180 s.
- **"Is the plastic zone clear of the boundary?" cannot be a fraction of the
  domain extent on a GRADED mesh** — the outermost hex here is 3.1 m wide with
  its centroid at 8.45 m of 10, so an element touching the roller passes a
  "< 90 %" test. Test element-COLUMN membership and report
  `n_yielded / n_in_column`.
- **`TaskStop` kills the bash wrapper, not the Python child** — but a
  foreground Bash *timeout* DOES kill the child. Check `Get-CimInstance
  Win32_Process` either way.

## Reproducing

```bash
C:/Users/nmb/AppData/Local/Programs/Python/Python312/python.exe Ladruno_files/testbed/hypo_bearing/bearing_backbone.py all
```

```bash
C:/Users/nmb/AppData/Local/Programs/Python/Python312/python.exe Ladruno_files/testbed/hypo_bearing/dp_collapse.py all
```

`dp_analyze.py` needs matplotlib, which only the apeGmsh venv has — it sets
`ADR79_NO_ENGINE=1` so it never loads the solver and the venv's hijacking
`.pth` is harmless.

```bash
C:/Users/nmb/venv/opensees_env/Scripts/python.exe Ladruno_files/testbed/hypo_bearing/dp_analyze.py
```

`.ladruno` recorder output is gitignored — `viewer_run.py` plus the committed
`bearing_mesh.npz` regenerate it exactly.
