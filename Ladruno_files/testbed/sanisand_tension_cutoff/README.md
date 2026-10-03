# WP-152 — tension cutoff (separation): the oracle-side definition

Plan and results: [`Ladruno_implementation/152_sanisand_tension_cutoff.md`](../../../Ladruno_implementation/152_sanisand_tension_cutoff.md).

- `tc_oracle.py` defines the NORMAL ⇄ SEPARATED state machine. It is a driver around the **unmodified** WP-134
  oracle: the WP-151 testbed copy `../sanisand_reseat_r1/sanisand_r1`, which carries the R1 toggles.
  - The oracle's `p_floor` stop, where the exact trajectory reaches p → 0 inside an increment, is the counterpart of
    SAS-ME's tension refusal (trigger E1).
  - The oracle integrates exactly, so it has no accuracy or cost failures. Trigger E2 (codes 4/9 at p0 < p_sep)
    exists on the C++ side only.
- **Runs** (CPython 3.11 with numpy + scipy):

  | command | what it does |
  |---|---|
  | `py -3.11 tc_oracle.py` | the three element paths (isotropic, triaxial extension, 4 open/close cycles) from an isotropic 2 kPa state: events and net work |
  | `py -3.11 tc_oracle.py --fixture <worktree>` | writes `tests/data/wp152_oracle_paths.json` from the post-flip state the C++ driver recorded (`tests/data/wp152_flip_state.json`, `wp152_cutoff_tools.record_flip_state()`) |

- The C++ gate is `tests/test_ladruno_sanisand_tension_cutoff.py`, with the helper `tests/wp152_cutoff_tools.py`.
