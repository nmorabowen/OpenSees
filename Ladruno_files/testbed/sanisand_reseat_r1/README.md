# WP-151 — R1 oracle prototype: the SANISAND α_in re-seat singularity

Memo: [`Ladruno_implementation/151_sanisand_reseat_singularity.md`](../../../Ladruno_implementation/151_sanisand_reseat_singularity.md).

- `sanisand_r1/` is a **copy** of WP-134's oracle (`Ladruno_scripts/sanisand_reference` @ 6cef73cc8). It adds
  the R1 mechanisms as `Options` toggles, **all OFF by default**. With them OFF it is bit-identical to the
  committed oracle (`check_identity.py`). The committed package is not edited.
- `r1_vs_wp134.diff` is the complete change.

| toggle (`Options`) | meaning |
|---|---|
| `h_reg="max"`, `h_eps=ε` | h = b0/max((α−α_in):n, ε), i.e. the C++ `-sasHFloor c_A` with ε = c_A·√(2/3)m |
| `h_reg="add"` | h = b0/(⟨a⟩ + ε), PM4Sand's C_γ1 linearised (comparison only) |
| `h_reg="max_soft"` | the floor only where b:n ≤ 0, i.e. WP-150's first form (comparison only) |
| `h_soft_kappa=κ` | where b:n < 0, h ≤ (1−κ)X/(⅔p\|b:n\|), i.e. `-sasSoftCap κ` |
| `reseat_delta=δ` | α_in re-seats only when (α−α_in):n < −δ, i.e. `-sasReseatHyst c_rev` with δ = c_rev·√(2/3)m |

Variant names in `r1common.variants()`: `B<c_A>` = floor, `T<c_rev>` = hysteresis, `S` = cap (κ 0.5),
`R150` = WP-150's form. For example, `T1B1S` = c_rev 1, c_A 1, κ 0.5 is the recommendation.

## Runs (CPython 3.11 with numpy + scipy; 6 worker processes)

| script | what | output |
|---|---|---|
| `check_identity.py <worktree>` | toggles OFF == the committed oracle (13 cases) | stdout |
| `extract_refusers.py` | one-off: the five wall refusers from the WP-138 Esmeralda checkpoints → `data/refuser_states.csv` (bit-identical rebuild) | `data/` |
| `a3_fan.py [ndir] [variants…]` | (a) the wall fan: 5 states × 32 directions × {3e-6, 3e-5} | `out/fan.json` (merged) |
| `a12_ring_reproducer.py [variants…]` | (a) b8 1950/3, 1950/2 and the WP-128 reproducer | `out/a12*.json` |
| `zeno_trace.py` | the Zeno re-seat sequence at E_B 1880/1 | `out/zeno_trace.json` |
| `b_tests.py`, `b_run.py [--variants a,b] [tests…]`, `b_analyse.py` | (b) calibrated behaviour: monotonic + cyclic element tests | `out/b_metrics.json`, `out/b_summary.md` (histories not committed, ~100 MB) |
| `cyc_pilot.py`, `cyc_sensitivity.py`, `cyc_gate.py`, `cyc_toyoura.py` | the CTXu gate: DM04's own perturbation sensitivity, every variant × perturbation, the c = 0.80 control, Toyoura | `out/cyc_*.json` |
| `fan_c080.py` | the wall fan with the Lode parameter c = 0.80 vs 0.71 (DM04), on the SAME c = 0.71 states. This measures c-sensitivity only; the c = 0.80 footing walls anyway (`fan_c080_bvp.py`) | `out/fan_c080.json` |
| `fan_fh_wall.py [--extract] [variants…]` | why floor + hysteresis WITHOUT the cap walls on the footing: fh's 11 NonPosH points (step 317), DM04 / fh / fhc and κ 0.25 / 0.75, with H/X under the floored h (states in `data/refuser_states_fh.csv`) | `out/fan_fh_wall.json` |
| `m1_small_strain.py` | review #893 M1: R1 vs DM04 at the SAME small strain (CVSS γ 1e-6…3e-4, undrained TC ε_a 1e-6…1e-4, 10 CVSS cycles at γa 1e-5 with c_rev 1 and 2) | `out/m1_small_strain.json` |
| `fan_c080_bvp.py [--extract] [variants…]` | the wall fan on the c = 0.80 footing's OWN wall states (`C080_EB_off`, step 308), at c = 0.80: DM04 vs R1 (states in `data/refuser_states_c080.csv`) | `out/fan_c080_bvp.json` |
| `c_jitter.py [c1 c2 c3]` | (c) jitter chains, objectivity (λ → 0), trial-direction continuity | `out/c1.json`, `c2.json`, `c3.json` |
| `r1plots.py [zeno fan mono gate continuity objectivity]` | the memo's figures | `out/fig/` |
| `cxx_fan.py <bin> <out> [flags…]` (CPython 3.12 `-S`) | the same fan on a C++ build | `out/cxx_fan_before.json` (pre-WP-151 binary) |
| `cxx_check_wp151.py <worktree>` (3.12 `-S`) | the WP-151 build: byte-identity, refusals per prototype, vs the oracle | `out/cxx_check_wp151.txt` |
| `make_byteid_baseline.py <bin>` (3.12 `-S`) | records `tests/data/wp151_sasme_byteid_baseline.json` on a build WITHOUT WP-151 | — |
| `cxx_fd_tangent.py <worktree> [bin] [wall ring cap kink]` (3.12 `-S`) | tangent consistency where R1 acts: Richardson FD of tiny replays vs the returned continuum tangent (floor, ring rows, cap-binding states, the cap's switch-on kink) | stdout |

`out/cyc_gate_v1_compounded.json` and `out/cyc_toyoura_v1_compounded.json` are the first gate runs. Their
perturbation wrapper compounded across tasks in one worker, so their per-run labels are approximate. They are
kept as the record of that correction. The clean reruns are `cyc_gate.json` and `cyc_toyoura.json`.

`R1_literature_survey.md`: the full literature survey, with sources and tags.

The C++ gate is `tests/test_ladruno_sanisand_reseat_r1.py` (helper: `tests/wp151_reseat_tools.py`).
