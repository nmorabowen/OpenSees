# External MPM reference: Kratos Multiphysics (WP-253)

One elastic, explicit 3D MPM case run in a permissively licensed code, to put the fork's mode-F
throughput beside it. This is a **reference point, not a validation**: the two codes do not
implement the same algorithm (see "Not matched").

## Install (isolated venv, outside the repo)

```
py -3.12 -m venv kratos_venv
kratos_venv\Scripts\python -m pip install KratosMPMApplication
kratos_venv\Scripts\python Ladruno_scripts\mpm_bench\external\kratos_column.py --N 10000 --threads 1
```

Versions used for the smoke: Python 3.12.10, `KratosMultiphysics` / `KratosMPMApplication` /
`KratosLinearSolversApplication` **10.4.4**, numpy 2.5.5 (Windows wheels, `cp312`). The MPM
package was `KratosParticleMechanicsApplication` up to v9; from v10 it is `KratosMPMApplication`
(`KratosMultiphysics.MPMApplication`). Pure `pip install KratosMPMApplication` pulls the core.

## Run

`kratos_column.py --N {1000,10000,100000} --threads n [--steps 20] [--warmup 5] [--scheme musl|usl|usf|central_difference] [--law ...] [--workdir dir]`
prints one JSON line (versions, N, dt, timing, a sanity block). N snaps to the fork's cell counts:
m = 9 gives N = 11,664 (the fork's "10^4"), m = 18 gives N = 93,312 ("10^5"), m = 4 gives 1,024.

**Timing window**: the wall of the full per-step cycle `AdvanceInTime` + `InitializeSolutionStep`
+ `Predict` (element search, mapping particle to grid) + `SolveSolutionStep` (grid solve, G2P,
stress update) + `FinalizeSolutionStep`, summed over the timed steps. It excludes mdpa writing,
reading, particle generation (`Initialize`, reported as `init_s`) and all output. Throughput =
N x steps / window. The fork's number is N x 20 / sum of the profiler `step` scope, which also
covers its whole step (support search through commit), so the windows are comparable in extent.
The `--threads` flag sets `OMP_NUM_THREADS` and `ParallelUtilities.SetNumThreads` (`kratos_threads`
in the JSON confirms it).

## Matched

| Item | Fork (`mpm_m0b/models.py`) | Kratos case |
|---|---|---|
| geometry | quarter column a x a x 2a, a = 1 m | same |
| background grid | cell h = a/m, 4m x 4m x (2m+2) cells | same hex8 grid, same cell count |
| particles | 8 per cell (2 per axis) | `MATERIAL_POINTS_PER_ELEMENT` = 8 on hex8 (2x2x2) |
| material | E = 1e6 Pa, nu = 0.3, rho = 2000 | same values |
| gravity | 9.81 down | `MP_VOLUME_ACCELERATION` / `MP_ACCELERATION` |
| symmetry | `fix xmin x`, `fix ymin y` | grid nodes fixed in X at x = 0 and Y at y = 0 |
| time step | 0.25 dt_crit, dt_crit = h / sqrt(M/rho), M = constrained modulus | same formula, set as `time_step` |
| integration | explicit, MUSL | `forward_euler` + `stress_update: musl` (explicit lumped mass) |
| I/O in window | none | none |

## Not matched (read before quoting a ratio)

1. **Shape functions.** The fork uses GIMP (support wider than a cell). Kratos MPM has linear
   hexahedral grid shape functions (standard MPM, support = the cell's 8 nodes). No GIMP or
   B-spline in Kratos. Kratos therefore does *less* transfer work per particle (8 nodes against
   up to 64), so the comparison is, if anything, favourable to Kratos.
2. **Constitutive model.** Fork: `ElasticIsotropic` driven by the particle strain increment.
   Kratos: `HyperElasticNeoHookean3DLaw` (finite strain; this wheel registers it, not
   `HyperElastic3DLaw` / `LinearElastic3DLaw`, which are classes but are not registered
   under those names). Same E and nu; the stress differs at finite strain, negligibly for the
   ~1e-4 per-step strains here. `--law LinearElasticIsotropic3DLaw` is accepted by the script.
3. **Floor.** The fork has a Coulomb wall (mu = 0.4) on zmin. The Kratos case fixes the vertical
   displacement of the z = 0 grid nodes (frictionless, no penetration). Its cost is a handful of
   nodes; the 25-step window barely moves the column.
4. **Time-step control.** Both use a fixed dt of 0.25 dt_crit computed from the same formula;
   Kratos does not adapt it.
5. **Threading.** Kratos threads with OpenMP inside its strategy; the fork's mode F is serial in
   the baseline. Thread scaling is not a like-for-like comparison before WP-253's speedup lands.
6. **Build.** Kratos is the official release wheel (MSVC, its own optimisation flags). The fork is
   our MSVC/MKL build. No attempt was made to rebuild either.
7. `Initialize` includes Kratos' bin-based search structure and MP generation; `Predict` runs the
   element search every step. The fork's support search is inside its step too. Not separated.

## Smoke result (one repeat, this machine, i7-13700KF, quiet but not isolated; NOT the pre-registered protocol)

| N | threads | window (20 steps) | median step | particle-steps/s |
|---|---|---|---|---|
| 11,664 | 1 | 2.293 s | 114.6 ms | 1.02e5 |
| 11,664 | 8 | 0.570 s | 29.3 ms | 4.09e5 |

Fork (232g M2, 1 thread, 3D elastic, 10^4): 1.368e6 particle-steps/s, i.e. about 13x Kratos at
1 thread. Sanity: mean particle vertical velocity at the end of the window is -0.219 m/s at
both thread counts (free fall for the 25-step 0.061 s window would be -0.60), max |vz| 0.263;
particle count unchanged (no particles lost).

## Licence and acknowledgement

Kratos Multiphysics is BSD licensed (PyPI metadata for `KratosMPMApplication` 10.4.4 reads
`License: BSD-4-Clause`; the wheel ships no licence file; the text is `license.txt` in
https://github.com/KratosMultiphysics/Kratos). Only the installed release wheels were run; no
Kratos source was copied into this repository. The script uses Kratos' public Python API and
writes its own mdpa files; nothing from it is to go into apeGmsh.

Acknowledgement: Kratos Multiphysics (CIMNE and the Kratos Team) and the authors of its
MPMApplication; see the application's README in the Kratos repository for the papers to cite.

## Published numbers (not matched, for context only)

Not collected: the Kratos run succeeded, so the fallback is not needed.
