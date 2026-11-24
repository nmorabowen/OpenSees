# mpm_bench — the WP-253 MPM performance harness

Implements sections 1-4 of the pre-registration
`Ladruno_implementation/253_mpm_performance.md` (issue #1097). It REUSES the M0b harness
(`Ladruno_scripts/mpm_m0b/common.py`, loaded as `m0b_common` by `benchlib.py`): a fresh child
process per run (`--child cfg.json`), thread env pinned before import, affinity mask (nmb:
`0xFFFF`, the P-cores), the profiler HDF5 reader, `median_spread`, the M1 mode-A decomposition.
Model constants (E, nu, rho, g, DP cone, block ramp) are read from `mpm_m0b/models.py`.

| file | what |
|---|---|
| `models.py` | the column (mode F) and block (mode A), 2D and 3D, three materials, two bases |
| `bench.py` | the matrix, child, parent loop; subset flags `--modes --bases --mats --dims --N --threads --repeats --list --dry-run --smoke --digest --no-deep --rerun-from` |
| `summarize.py` | median / spread, the >10 % re-run flag, "within noise", S(n), G(n), markdown tables |
| `regress.py` | the five-row regression check against `baseline_nmb.json` |
| `esmeralda/wp253_bench.sbatch` | the 10^6 and 24-thread rows (not submitted by anything here) |

## Running

```
python3.12 Ladruno_scripts/mpm_bench/bench.py --dry-run                  # the 81 nmb rows, model sizes
python3.12 Ladruno_scripts/mpm_bench/bench.py --smoke --N 1000,10000     # smoke, 1 repeat
python3.12 Ladruno_scripts/mpm_bench/bench.py --modes F --dims 3 --N 100000 --threads 1,8
python3.12 Ladruno_scripts/mpm_bench/summarize.py <out>/bench_results.json [--baseline base.json] --md t.md
python3.12 Ladruno_scripts/mpm_bench/bench.py --rerun-from <out>/bench_results.json   # 5-repeat rerun of spread > 10 % rows
python3.12 Ladruno_scripts/mpm_bench/regress.py [--record]
```
`--pyd` / `--exe` default to `dist/bin` of this checkout. Results go to `results/` (per-run files
git-ignored; `*_results.json` aggregates may be committed). `--digest` records a sha256 per particle
array after the timed steps; `summarize` compares them across thread counts (exact equality).
`--site esmeralda` selects the Linux matrix (mode F, N up to 10^6, threads 1, 4, 8, 24).
No build or other load may run during a timed run.

## Provenance (every row)
`ladrunoBuild` of `OpenSees.exe` (run on a one-line Tcl script) and of the pyd (a probe child);
**the invocation is refused if they differ** (and a row whose child reports another build is
marked failed). Also machine (CPU model, logical cores, RAM, OS), MKL (file version of the
`mkl_core`/`mkl_rt` the process loaded, `mkl_get_version_string` when `mkl_rt` is mapped, else the
DLLs beside the binary; this build has no `mkl_rt`, so the file version, e.g. 2026.1.0.0, is
what is recorded), `ladrunoThreads` read back, `OMP/MKL_NUM_THREADS`, N, steps, Newton k (mode A).

## Model choices and why
- **Mode F column.** 3D: quarter column, body `a x b x 2a` (a = 1), `fix xmin x`, `fix ymin y`,
  `mpm wall zmin 0.4`; 2D plane strain: half column `a x 2a`, `fix xmin x`, `mpm wall ymin 0.4`
  (the vertical axis is y in 2D). Grid 4x wider than the column so the collapse never leaves it.
  ppc 2 per axis, dt = 0.25 dt_crit, as M2. 5 warm-up + 20 timed steps in one `mpm analyze`.
- **N within 10 %.** M2's cube quarter column (N = 16 m^3) cannot hit 10^4 within 10 %
  (8192 / 11664). So the 3D column keeps height = 2 x width in x and lets the y extent
  `b = ny h` differ from a (`models.column_cells`): N = 16 nx^2 ny (10368, 98496, 998400 at the
  three sizes). 2D: N = 8 nx^2 (9800, 100352, 1002528). Block: 8 n^3 (10648, 97336, 10^6) and
  4 n^2. The actual N is recorded.
- **Mode A block.** M1's: confined block, gravity ramp to 0.75 over 8 steps, PARDISO,
  `emptyNodes renormalize`, `massCutoff 5e-2`, NormDispIncr 1e-6; 2D: `ymin xy`, `xmin/xmax x`.
  2 warm-up + 6 timed, one profiler run per step. T uses the `mpm.a.step` wall; also
  T per Newton iteration (N k / W) and M1's t_rb / t_it / R / f.
- **Materials.** ElasticIsotropic and DruckerPrager in SI as M0b (DP: column phi 30 deg, sigmaY 10 Pa;
  block rho 0.21, sigmaY 400 Pa). LadrunoSANISAND: the fork's kPa deck
  (`tests/test_ladruno_mpm_soil_admission.py`), `IntScheme 2`, rho 2.0 t/m3, g 9.81, same geometry;
  mode F only on bspline2 (232k K0), mode A 3D gimp 10^4 1 thread. It starts in stage 0
  (elastic) for the warm-up, then `mpm stage 1 1` before the timed steps; dt is lowered after the flip if
  the critical step dropped (never above 0.25 of it). Its plastic call dominates (~90 % of a step), so the
  SANISAND rows say little about the transfer loops; the pre-registration says so too.
- **Phase split.** Shares of the `step` wall for the scopes listed in section 2, found anywhere
  below `step`; `other` = remainder. A build without them reports `mpm.advance`, `commit`, `other`
  (`scopes_present: false`). The material share comes from ONE separate deep run per row
  (`profile -deep`, `mpm.f.material` elem_by_type, M2's method); `--no-deep` skips it.
- **Noise.** 3 repeats, median, spread = (max - min) / median; > 10 % with < 5 repeats is flagged
  (`*`), `--rerun-from` re-runs those with 5 repeats and `summarize` takes the file with the most
  repeats per row. Two medians differ only if |ratio - 1| > sum of spreads.
- **regress.py.** Rows: 3D gimp elastic 10^5 F at 1 and 8 threads; 3D gimp DP 10^5 F at 1; 2D
  bspline2 elastic 10^5 F at 8; 3D gimp elastic 10^4 A at 1. 3 repeats, no deep run. A row
  regresses when median T falls by more than max(10 %, 2 x baseline spread). Exit 0 / 1 / 2.
  The orchestrator records `baseline_nmb.json` with `--record` (not committed here).
- **esmeralda.** `sbatch --exclude=node1,node2,node12,node13 --mem-per-cpu=4G --exclusive ...`,
  `OMP_PROC_BIND=close OMP_PLACES=cores`, runs from a fresh clone given as an argument and refuses
  `~/ladruno_build_test/OpenSees`; `--affinity none`; machine and `lscpu` saved with the results.
  Mode A is nmb-only (as the pre-registration's matrix).
