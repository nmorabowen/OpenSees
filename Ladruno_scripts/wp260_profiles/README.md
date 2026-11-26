# WP-260 Wave 0 (`260-w0`) — Phase-0 Profiler Baseline

- **Binary Build SHA (`ops.ladrunoBuild()`):** `7e91aaf1c48870baaa9720fbbd6d5ceca669da72` (pre-`#1147`/`#1153`/`#1154` baseline)
- **Host:** i7-13700KF (`8` P-cores / `16` logical HT pinned via affinity `0xFFFF`, zero concurrent builds)
- **Protocol:** `ops.profiler('start', '-deep', '-perStep')` at 1T (`OMP_NUM_THREADS=MKL_NUM_THREADS=1`, `ladrunoThreads(1)`) + coarse `ops.profiler('start', '-perStep')` at 8T (`OMP_NUM_THREADS=MKL_NUM_THREADS=8`, `ladrunoThreads(8)`), `3` interleaved reps.
- **Accounting rule:** Every row sums disjoint scopes to 100.0% of profiled step time; `elem.tangent` is strictly split into element kernels (`elem_by_type` inside `elem_kern`) and sparse scatter (`soe.addA`).

## 1. Measured Phase-Share Table (1T `-deep`) & Production Wall Time (8T)

| Deck | Family | Steps | `soe.sym` | `soe.fac` | `soe.tri` | `soe.addA` | `elem_kern` | `newStp_o` | `unbal_o` | `upd_o` | `dc_o` | `cmt_el` | `cmt_rec` | `other` | 1T med [min..max] | 8T med [min..max] | 1T/8T |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `d1a_brick_sanisand` | D1a-Brick-SANISAND | 24 | 0.0% | 0.1% | 0.0% | 0.4% | 98.4% | 0.5% | 0.1% | 0.1% | 0.0% | 0.2% | 0.0% | 0.1% | 47.026s [40.210..61.353] | 39.572s [39.312..60.832] | 1.19x |
| `d1b_brick_concrete3d` | D1b-Brick-Concrete3D | 16 | 0.1% | 1.7% | 0.8% | 0.5% | 61.2% | 6.4% | 0.2% | 28.9% | 0.2% | 0.0% | 0.0% | 0.1% | 10.199s [10.102..14.263] | 10.632s [10.236..13.863] | 0.96x |
| `d2_cpml_wave` | D2-CPML-Implicit | 30 | 2.3% | 26.8% | 4.1% | 6.3% | 45.2% | 6.3% | 1.4% | 1.2% | 0.6% | 4.0% | 0.0% | 2.0% | 5.224s [5.086..5.441] | 4.680s [4.547..4.917] | 1.12x |
| `d3a_explicit_cdl_cpml` | D3-CentralDifferenceLadruno | 100 | 0.0% | 0.0% | 0.0% | 0.0% | 57.6% | 19.7% | 1.1% | 1.0% | 0.5% | 19.6% | 0.0% | 0.4% | 8.327s [7.856..9.019] | 7.809s [7.704..9.139] | 1.07x |
| `d3b_explicit_bathe_cpml` | D3-ExplicitBathe | 80 | 0.0% | 0.0% | 0.0% | 0.0% | 52.8% | 18.8% | 0.4% | 19.4% | 0.1% | 7.9% | 0.0% | 0.6% | 17.987s [14.362..18.815] | 17.954s [14.101..19.119] | 1.00x |
| `d3c_explicit_sms_cons` | D3c-SMS-Consistent | 250 | 3.4% | 0.4% | 16.8% | 0.1% | 64.1% | 4.6% | 2.2% | 3.8% | 2.1% | 1.8% | 0.0% | 0.9% | 3.264s [3.088..3.656] | 3.146s [2.658..4.012] | 1.04x |
| `d4_mpm_modeA` | D4-MPM-ModeA | 12 | 8.8% | 8.8% | 1.7% | 9.7% | 38.1% | 0.6% | 3.3% | 2.3% | 19.4% | 1.6% | 0.0% | 5.6% | 1.975s [1.918..2.019] | 1.677s [1.662..1.768] | 1.18x |
| `d5_mcft_2cld12` | D5-MCFT-2CLD12 | 24 | 0.0% | 0.0% | 0.1% | 0.0% | 97.8% | 0.0% | 0.0% | 0.0% | 0.0% | 2.0% | 0.0% | 0.0% | 9.728s [7.467..10.516] | 9.816s [7.778..10.235] | 0.99x |

### Element Kernel Composition (`elem_by_type`) & SANISAND Substep Share
- **`d1a_brick_sanisand`**: LadrunoBrick (33002)=98.4% *(SANISAND `sasME.substeps`: 87.7%)*
- **`d1b_brick_concrete3d`**: LadrunoBrick (33002)=61.2%
- **`d2_cpml_wave`**: LadrunoCPML (33023)=32.7%, stdBrick (56)=12.5%
- **`d3a_explicit_cdl_cpml`**: LadrunoCPML (33023)=54.0%, LadrunoBrick (33002)=3.7%
- **`d3b_explicit_bathe_cpml`**: LadrunoCPML (33023)=49.8%, LadrunoBrick (33002)=3.0%
- **`d3c_explicit_sms_cons`**: LadrunoBrick (33002)=64.1%
- **`d4_mpm_modeA`**: LadrunoMPMPoint (33026)=38.1%
- **`d5_mcft_2cld12`**: LadrunoForceBeamColumn3d (33025)=97.8%

## 2. Phase -> WP-260 Package Map

| Scope / Phase | Mapped Package(s) | Status on `7e91aaf1c` |
|---|---|---|
| `soe.factor`, `soe.symbolic`, `soe.trisolve`, `soe.addA` | `260a` (PARDISO/MUMPS hygiene), `260d` (reuse), `260g` (SP/IR), `260m` (AMG) | `260a` / `260d` / `260g` unblocked |
| `LadrunoCPML` (`33023`) + explicit `update`/`residual`/`commit` (`D3a`/`D3b`/`D3c`) | `260b` (Explicit & Wave Hygiene) | Unblocked |
| `elem_by_type` (`LadrunoBrick` SANISAND/Concrete3D, `LadrunoForceBeamColumn3d` MCFT, `LadrunoMPMPoint`) | `260f` (Threading+FastScatter), `260c2` (Concrete3D stress-only/implex), `260h` (Analytic return maps), `260k` (SoA) | `260f` blocked (`#1147/#1153`); `260c2`/`260h` blocked (`#1153/#1154`) |
| Cross-cutting OMP / Allocator / LTO (`1T/8T <= 1.0x` on threaded/MKL decks) | `260e` (D4: `libiomp5md` unification + `mimalloc` + `/GL`) | Unblocked (benchmark after Wave-1 isolation) |
| `commit.record` | Outside WP-260 (recorders; flagged if >=10%) | Observed `0.0%` on all D1 baseline decks (no heavy file recorders enabled in bench decks) |

## 3. Plan §0.5 Amdahl Ranking (`p >= 10%` and `G = 1/((1-p)+p/s) >= 1.10x`)

| Rank | Pkg | Tier | Status / Blocker | Target Scope | Conservative `s` | Peak Deck | Peak `p` | Peak `G` | Gate |
|---:|---|---|---|---|---|---:|---:|---:|:---:|
| 1 | `260f` | G-BYTE/G-EQUIV | Blocked (#1147/#1153) | elem_kernel + soe.addA across D1a/D1b/D4/D5 (currently serial on 8 P-cores) | 2.50x | `d1a_brick_sanisand` | 98.8% | **2.46x** | PASS |
| 2 | `260l` | G-OPT | Post-Wave 1 | Whole explicit step (D3a/D3b) | 2.00x | `d3a_explicit_cdl_cpml` | 100.0% | **2.00x** | PASS |
| 3 | `260h` | G-EQUIV/G-OPT | Blocked (#1153/#1154) | elem_kernel on D1a SANISAND / D1b Concrete3D | 2.00x | `d1a_brick_sanisand` | 98.4% | **1.97x** | PASS |
| 4 | `260k` | G-OPT | Post-Wave 1 | elem_kernel on D1a SANISAND | 2.00x | `d1a_brick_sanisand` | 98.4% | **1.97x** | PASS |
| 5 | `260b` | G-BYTE/G-EQUIV | Unblocked | LadrunoCPML + explicit update/residual/commit (D3a/D3b/D3c) | 1.80x | `d3b_explicit_bathe_cpml` | 96.4% | **1.75x** | PASS |
| 6 | `260c2` | G-BYTE/G-EQUIV | Blocked (#1153/#1154) | LadrunoBrick (33002) on D1b Concrete3D | 2.00x | `d1b_brick_concrete3d` | 61.2% | **1.44x** | PASS |
| 7 | `260a` | G-BYTE+G-OPT | Unblocked | soe.symbolic + soe.factor + soe.trisolve + soe.addA | 1.70x | `d2_cpml_wave` | 39.5% | **1.19x** | PASS |
| 8 | `260d` | G-EQUIV/G-OPT | Unblocked | soe.factor + elem_tangent_kernel on implicit transient (D1a/D2/D4) | 1.60x | `d2_cpml_wave` | 37.6% | **1.16x** | PASS |
| 9 | `260e` | G-BYTE | Unblocked (cross-cutting build bench) | Global build/runtime allocator + single OMP pool across all decks | 1.15x | `d1a_brick_sanisand` | 100.0% | **1.15x** | PASS |
| 10 | `260m` | G-OPT | Post-Wave 1 | soe.factor + soe.trisolve + soe.symbolic on 3D implicit (D1a/D1b/D2) | 1.50x | `d2_cpml_wave` | 33.2% | **1.12x** | PASS |
| 11 | `260g` | G-OPT | Unblocked | soe.factor + soe.trisolve on FP64-safe implicit (D1a/D1b/D2/D4) | 1.50x | `d2_cpml_wave` | 31.0% | **1.12x** | PASS |
| 12 | `260c1` | G-BYTE/G-EQUIV | Unblocked | LogStrain / UP (finite-strain suite; N/A on small-strain D1) | 3.00x | `-` | 0.0% | **1.00x** | BELOW |

## 4. Proposed Wave-1 Picks (for Owner Confirmation & Claude §4 Review)

### Top Overall by Projected `G` (shows blocked dependencies)
- **`260f`** (G-BYTE/G-EQUIV, **Blocked (#1147/#1153)**) — `3D/MPM/MCFT Threading + Persistent Assembly + LadrunoFastScatter`: peak `p = 98.8%` on `d1a_brick_sanisand`, conservative `G = 2.46x` at `s = 2.50x`.
- **`260l`** (G-OPT, **Post-Wave 1**) — `LadrunoExplicitEngine Flat Explicit + CPML + Clustered LTS`: peak `p = 100.0%` on `d3a_explicit_cdl_cpml`, conservative `G = 2.00x` at `s = 2.00x`.
- **`260h`** (G-EQUIV/G-OPT, **Blocked (#1153/#1154)**) — `Analytic Return Maps (SANISAND SAS-RK, Concrete3D Cardano, Anderson/LBFGS)`: peak `p = 98.4%` on `d1a_brick_sanisand`, conservative `G = 1.97x` at `s = 2.00x`.
- **`260k`** (G-OPT, **Post-Wave 1**) — `LadrunoGPStore SoA + Active-Set Compaction + 7x7 Schur SANISAND CPPM`: peak `p = 98.4%` on `d1a_brick_sanisand`, conservative `G = 1.97x` at `s = 2.00x`.

### Recommended Immediate Unblocked Wave-1 Launches (2–3 Parallel WPs per §0.5)
- **`260b`** (G-BYTE/G-EQUIV, Unblocked): peak `p = 96.4%` on `d3b_explicit_bathe_cpml`, projected `G = 1.75x`.
- **`260a`** (G-BYTE+G-OPT, Unblocked): peak `p = 39.5%` on `d2_cpml_wave`, projected `G = 1.19x`.
- **`260d`** (G-EQUIV/G-OPT, Unblocked): peak `p = 37.6%` on `d2_cpml_wave`, projected `G = 1.16x`.

> **Dependency note:** Higher-ranked element/constitutive WPs (`260f`, `260l`, `260h`) wait for open upstream/fork PRs `#1147` / `#1153` / `#1154` (which add MPM & SANISAND/Concrete3D threading) before cutting their WP branches so `260-w0` re-ranks cleanly after those merges.
