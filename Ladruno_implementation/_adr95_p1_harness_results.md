# ADR-95 P1 harness — VERDICT

**Status:** harness ready, tolerant of the pre-token binary. Not a physics result.

## What was added (`Ladruno_files/testbed/hypo_bearing/quad_path_diag.py`)
- `--branch`: at every sampling station (the existing `--cond`/`--cond-at` cadence,
  UNIONED with `--dense LO HI STEP`, default `0.0108 0.0113 5e-6` = note 82 §7.3's
  bracket), loops every element/GP, pulls `ops.eleResponse(e,'material',gp,
  'ladrunoBranch')`. Records histogram (n0..n3), n_forced, detAmin/I1 extrema,
  gamma0/1 max, the 5 lowest-detAmin GPs and first 5 corner-branch GPs (element
  centroid as an approximate GP location) to `qpd_<tag>_branch.csv` (documented
  header comment; carries s/B, q, sigma_min/cond alongside the branch columns)
  and full per-GP arrays per station to `qpd_<tag>_branch.npz`. An always-fired
  final station is taken AT THE WALL (mirrors the existing diag-3 idiom).
- `--forensics`: on the first failed ladder attempt, and again at the step that
  hits the floor, records every ladder rung's `ops.testNorms()`/`testIter()`
  and the 10 largest `|nodeUnbalance|` DOFs with node coordinates (read after
  the failed `analyze(1)`, i.e. at the reverted last-converged state) to
  `qpd_<tag>_forensics.json`.
- Empty/short `ladrunoBranch` responses are counted, never raised — verified
  live (below): the harness runs unchanged against a binary with no token.
- `h20_prandtl.py`: `_BIN` now honours env `ADR95_DIST`, falling back to the
  unchanged worktree-relative `dist/bin` when unset.

## GP counts used (`GP_COUNTS`, source-verified, `SRC/element/ladrunoBrick`)
h8bbar = 8, h8std = 8, h20uri = 8 (`NGPU`, uniform reduced 2x2x2),
h20std = 27 (`NGP`, full 3x3x3).

## Smoke test (against the MAIN CHECKOUT binary, `ADR95_DIST` override)
```
ADR95_DIST=C:/Users/nmb/Documents/Github/OpenSees/dist/bin \
  py -3.12 quad_path_diag.py --elem h8bbar --h0 1.0 --branch --forensics \
  --suffix _smoke --sfrac 0.002
```
`ladrunoBuild()` printed `0e5163e430...`. Ran to TARGET in ~2 s (13 steps, 0
failures): `qpd_..._branch.csv/.npz` and `_forensics.json` all written; 1600
`ladrunoBranch` calls, all empty, all counted (`0/1600` live, no crash).

Second run exercised the harder paths (`h20uri --cond --branch --forensics
--dense 0.0005 0.003 0.0002`, capped `--tmax 90 --budget 40`): 12 failed
attempts, `[forensics] captured 'first_fail'` with 3 ladder-rung norm
histories + top-10 node-unbalance table; 14 branch stations written (dense +
final); the AT-THE-WALL branch row carried real `sigma_min`/`cond` from the
existing diag-3 sampler (`3.65e-4` / `4.68e4`). Both smoke output sets deleted
after inspection (dry runs, not campaign data).

## P1 legs (plan §3, to run once the P0 token lands)
```
py -3.12 quad_path_diag.py --elem h20uri --h0 1.0 --branch
py -3.12 quad_path_diag.py --elem h8bbar --h0 1.0 --branch   # linear control
```
