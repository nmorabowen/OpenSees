# sanisand_lowp_clamp — WP-230 (TIMs F25)

The oracle side of `-sasLowP clamp`. Plan and results: `Ladruno_implementation/230_sanisand_lowp_clamp.md`.

| file | what |
|---|---|
| `clamp_oracle.py` | the clamp state machine as a definition, a driver around the unmodified WP-134 oracle (the R1 testbed copy `../sanisand_reseat_r1/sanisand_r1`). `mode="separation"` is WP-152's `../sanisand_tension_cutoff/tc_oracle.py`, record for record (`--selfcheck`). `--fixture` writes `tests/data/wp230_oracle_paths.json`: five element paths × the two α rules. Python 3.11 with scipy. |

The C++ side is `tests/test_ladruno_sanisand_lowp_clamp.py` (`tests/wp230_clamp_tools.py` holds the
byte-identity recorder).
