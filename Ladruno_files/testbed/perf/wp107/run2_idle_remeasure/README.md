# run2 — the idle-box re-measure (WP-107 red-team S4)

The raw data behind `../RESULTS.md` §2. It **supersedes** `../run1/curve_el_elastic_*.csv`
for the speed-up table; run1 is kept because it is what the first version of the PR
reported and the discrepancy is the point.

Command (once per `threads` ∈ {1,2,4,8} × `rep` ∈ {0,1,2}, from the worktree root):

```
set PYTHONPATH=<worktree>\dist\bin
set MKL_NUM_THREADS=1
python3.12 -u Ladruno_files/testbed/perf/wp107/wp107_strip_bench.py ^
    --threads <T> --h 0.05 --steps 12 --mat elastic --system Pardiso ^
    --out Ladruno_files/testbed/perf/wp107/run2_idle_remeasure/el_t<T>_r<R>.csv
```

**Box idle** — that is the whole difference from run1, which was taken while two sibling
wp/106 SANISAND jobs were running and was captioned "speed-ups are therefore lower
bounds". The caveat had the wrong sign: contention was suppressing the 8-thread
oversubscription penalty as much as it was inflating the serial baseline.

| threads | per-step wall (s, min of 3) | mean of 3 | speed-up | run1 reported |
|---|---|---|---|---|
| 1 | 0.16014 | 0.17301 | 1.00x | 0.16557 / 1.00x |
| 2 | 0.15550 | 0.16865 | 1.03x | 0.16518 / 1.00x |
| 4 | 0.15474 | 0.17034 | 1.03x | 0.15146 / 1.09x |
| 8 | 0.17050 | 0.19819 | **0.94x — a regression** | 0.14872 / **1.11x** |

run1's 1.11x is above this deck's own Amdahl ceiling of 1.08x (loop A = 7.80 % of step),
which is what should have flagged it as noise at the time. The red team measured
1.03 / 1.05 / 0.98x independently on the same deck; this run gives 1.03 / 1.03 / 0.94x.

## The `fieldmd5` column is new (red-team S6)

The harness used to gate bit-identity on `settlement_m` — a **prescribed** `sp` value,
identical by construction at every thread count — and `footing_load`, a **sum** of
reactions in which a per-element scramble can cancel. Interior displacements were never
compared. `fieldmd5` is an md5 **per step** over every node's `nodeDisp` and
`nodeReaction` at `repr()` precision plus one element's stress vector.

All 12 runs here share one rollup digest:

```
FIELDMD5 ... field=2ddbd45b00ffcc207b374600bbe43195
```

i.e. bit-identity holds across 1/2/4/8 threads on the full field, not just on two summary
numbers.
