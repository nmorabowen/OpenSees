# ADR-79 bearing campaign — session handoff (2026-07-30)

Status: **the campaign is complete and merged.** PRs #684, #686, #689 (this
line) plus #687/#688 (the concurrent locking line) are all on `ladruno`.
Nothing is in flight. This note exists so a new session can pick up the
*follow-on* work without re-reading the campaign.

## Where everything lives

| what | where |
|---|---|
| runner, probes, figure scripts, CSVs | `Ladruno_files/testbed/hypo_bearing/` |
| fork-side write-up | that folder's `RESULTS.md` and `README.md` |
| fork decision record | ADR 79 §8 (`79_ladruno_geom_hypo_adr.md`) |
| **project-facing analysis** | TIMs vault notes **16** and **17**, `C:\Users\nmb\Dropbox\obsidian\Ladruño\TIMs\Reference Model\` |
| traps | `LEDGER_quirks.md` (2 entries from this line) |

Worktree `C:\Users\nmb\Documents\Github\OpenSees-hypo`, branch
`claude/adr79-bearing-deeper`, built `dist/bin` current.

## The answer, in three lines

1. **The geometry ladder is exonerated.** −5.4 % (corot) to +7.5 % (hypo)
   against a `-geom linear` base, no limit point out to s/B = 15 %, tangent flat
   at 6847 kPa/m. ADR-61 **OQ8 answered negatively** — the finite-strain soil law
   it contemplated is not needed. OQ6's second item (`-geom corot` on
   `LadrunoUP`) is delivered; its first (relax the contact `ndf == 3` guard) is
   now the blocker for testing the **bonded interface**, the last kinematic
   candidate.
2. **The capacity anchor is wrong twice.** PDMY's failure surface is a measured
   Drucker–Prager cone (α constant to 3.9 % across the Lode extremes while the
   MC angle swings 31.5° → 54°), calibrated in triaxial compression; Vesić's
   `N_γ` is Mohr–Coulomb, and the plane-strain equivalent is 53.7° ⇒ 10.8 MPa
   against the 207 kPa anchor, a factor of 52. And the failure **mode** is
   punching, not the general shear `N_γ` describes.
3. **So the "9.3× over-strength" was never a well-posed anomaly.** The model
   sits at ~31 % of its own collapse load. The defect measurements all stand and
   still matter for *stiffness* — locking 34.3 % (#688), rate artifact 11.7 %.

## What to do next, in priority order

1. **Numerical collapse load with the actual cone.** This is the one calculation
   that converts §2 above from a strong hypothesis into a measurement. The 52×
   rests on Chen–Han plane-strain matching, which assumes *associated* flow that
   this non-dilatant set violates, and treats a square footing as plane strain.
   Cheap, and it gates the interpretation of every capacity ratio in the vault.
2. **Re-anchor or annotate the vault ratios** so nobody reads 9.3× as
   over-strength. Note 17 says this explicitly; the ratios themselves are fine
   as recorded data.
3. **A modelling decision for the project:** is a Lode-independent cone the
   intended constitutive assumption? It over-predicts plane-strain strength by
   ~50× on bearing capacity. Either calibrate φ for the regime of interest or
   use a Lode-dependent surface — a deliberate choice, not a solver matter.
4. **The interface**, once the contact `ndf == 3` guard is relaxed. It is the
   last untested kinematic candidate and this benchmark cannot speak to it (the
   footing is a smooth driven node patch, so an interface is absent by
   construction).
5. **A refined, unlocked mesh** — the only experiment that separates "this soil
   genuinely punches" from "four elements across B cannot form a band". Note
   that `-geom hypo` still refuses `-formulation bbar` (ADR 79 P2 reserved), so
   locking and large strain cannot yet be measured together on the coupled
   element; `corot` composes with `bbar` today.

## Traps this line paid for — do not rediscover

- **An installed Ladruno hijacks `import opensees`** via a site `.pth` that
  imports at interpreter startup, so a worktree `sys.path.insert` is a no-op and
  scripts silently run a build that refuses `-geom corot|hypo`. Use base Python
  3.12 and **assert `opensees.__file__`** (the runner does).
- **`NormDispIncr` is not evidence of equilibrium.** It reported convergence at
  σ_zz = +0.26 kPa *tension* under a 100 kPa compressive load once the tangent
  degenerated. Use `NormUnbalance` where that is possible.
- **A converged push increment is a property of the mesh**, not the problem.
  P3's 2.5 mm fails here for *every* geometry method including `-geom linear`;
  it is a transient first-loading shock off the stage-1 PDMY state. Drive with a
  2-point ramp + adaptive increment; KrylovNewton is the primary algorithm.
- `remove('sp', node, dof)` does **not** match a constraint made by `fix`.
  UmfPack `setSize` fails at 4 DOFs (use `FullGeneral`). `sinφ_mob` needs a
  **tension** guard as well as a confinement one. Field plots must be referenced
  to the gravity frame. Slab filters need a straddle test.
- **`Results.from_ladruno` is self-sufficient** — no `model.h5`, no deck↔FEMData
  tag reconciliation. `from_mpco` would have required both.
- **`TaskStop` kills the bash wrapper, not the Python child.** Three
  `field_run` processes once wrote the same `.npz` concurrently.

## Reproducing

```bash
C:/Users/nmb/venv/opensees_env/Scripts/python.exe Ladruno_files/testbed/hypo_bearing/build_mesh.py
```

```bash
C:/Users/nmb/AppData/Local/Programs/Python/Python312/python.exe Ladruno_files/testbed/hypo_bearing/bearing_backbone.py all
```

`.ladruno` recorder output is gitignored — `viewer_run.py` plus the committed
`bearing_mesh.npz` regenerate it exactly.
