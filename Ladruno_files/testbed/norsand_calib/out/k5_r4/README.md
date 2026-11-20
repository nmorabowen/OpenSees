# K5 strip legs on the round-4 engine (WP-144)

Engine `0d38e57b4` built on the Esmeralda login host. Each leg folder holds `steps.csv` (one row per committed step:
s/B, q [kPa], failed attempts, wall) and `summary.json`. The run logs (`log.log`) stay out of git (23 MB); they remain on
Esmeralda under `~/wp144/deck/runs/<leg>/logs/`.

`deck/` is the deck as run: `footing_ab_ns_k5.py` (a copy of `footing_ab_ns.py` with the K5 switches), its launcher
`launch_k5.py` and `leg_k5.sbatch`, one submit script per series, the comparators (`k5_*cmp.py`, `k5sum.py`), the
patches that made the copy (`k5_patch.py`, `patch_elem.py`, `patch_elem_ns.py`) and the NorSand parameter files.

Switches of the deck copy, read from the environment:

| Variable | Effect |
|---|---|
| `DP_ASSOC=1` | Drucker-Prager associated (psi = phi), cone matched to plane-strain Mohr-Coulomb |
| `FOOT_SMOOTH=1` | footing tied in y only (smooth) |
| `ELEM=<tech>` | `qbbar` (default), `qstd`, `qssp`, `qeas` (LadrunoQuad), `cst` (LadrunoCST), `lst` (LadrunoLST), `bez`, `bezb` (BezierTri6, std and `-bbar`); triangles 2 per cell, alternating diagonals |
| `ELEM_LST=1` | old name of `ELEM=lst` |

Series (prefix, submit script):

| Prefix | Content | Script |
|---|---|---|
| `K5r4_*` | NorSand, Kimura 1985 and Siddiquee 1999, three meshes, regular and sheared | `k5_rest_submit.sh` |
| `K5san_*` | SANISAND comparator at B/8 | `k5_rest_submit.sh` |
| `K5el_kimura_*` | NorSand with the elastic moduli scaled (G0x) | `k5_submit.sh` |
| `K5pi_*` | NorSand with pi_i0 of 1.0 to 2.0 p | `k5_pi0_submit.sh` |
| `K5sp_*` | NorSand, Siddiquee with the test pit as the box (14 B by 8 B) | `k5_sidpit_submit.sh` |
| `K5dp_*`, `K5dpa_*`, `K5dps_*` | Drucker-Prager phi 38: psi 0, associated, smooth footing | `k5_dp_submit.sh`, `k5_dp2_submit.sh` |
| `K5lsa_*`, `K5lsn_*` | Drucker-Prager on LadrunoLST: associated, psi 0 | `k5_lst_submit.sh` |
| `K5el_<tech>_b4/b8` | associated Drucker-Prager, eight element technologies | `k5_el_submit.sh` |
| `K5nse_*` | NorSand K3 on EAS quads and LadrunoLST | `k5_nsel_submit.sh` |

Readings are in the WP-144 notes (`state/00 — Current state.md`), not here.
