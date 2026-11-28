# WP-260 Wave 0 (`260-w0`) — SMOKE / PRE-#1153/#1154 VERIFICATION

> **Partial (D9): nothing can be picked from this table.** Present: none. Missing: `d1a_brick_sanisand[S]`, `d1a_brick_sanisand[L]`, `d1b_brick_concrete3d[S]`, `d1b_brick_concrete3d[L]`, `d2_cpml_wave[S]`, `d2_cpml_wave[L]`, `d3a_explicit_cdl_cpml[S]`, `d3a_explicit_cdl_cpml[L]`, `d3b_explicit_bathe_cpml[S]`, `d3b_explicit_bathe_cpml[L]`, `d3c_explicit_sms_cons[S]`, `d4_mpm_modeA[S]`, `d5_mcft_2cld12[S]`.

- **Binary SHA (verified `OpenSees.exe` == all `opensees.pyd` records):** `7e91aaf1c48870baaa9720fbbd6d5ceca669da72`
- **Environment (`OMP_*`/`KMP_*`/`MKL_*`/`LADRUNO_*`):** `{"LADRUNO_OPENSEESPY_ALIAS": "1", "LADRUNO_OPENSEES_QUIET": "1"}`
- **Affinity arg:** `none` | **Production threads:** `8`
- **Host CPU load gate (`--max-load`):** `<= 10.0%` (fail-closed on `cpu_load_before_pct`, `cpu_load_after_pct`, and `cpu_load_during_ext_pct` when present)
- **Per-run CPU load probes:** before: `min=13.3%, med=19.9%, max=38.7%` | after: `min=16.5%, med=23.4%, max=51.7%`
- **Step-1 exclusion (Finding 2):** step 1 runs BEFORE `profiler reset/start`, so step-1 symbolic/first-touch/D2 one-off factor costs are excluded from BOTH wall times AND rollup phase shares (and reported separately below).
- **Accounting & rankability gates (Findings 1 & 3, brief §6):** `residual = step_wall - accounted` is NEVER clamped (`sum > 100.5%` fails hard). Any deck with `deep_ovh < 0.95`, mode call-count mismatch, `|residual| > 0.5%`, unresolved remainder `> 5.0%`, or host CPU load `> 10.0%` (or missing) is marked `NO` (not rankable) and excluded from ranking.

## 1. Phase-Share Table (`1t_deep`) & 3-Mode Wall Times

| Deck[Size] | DOF~ | Steps | `soe.sym` | `soe.fac` | `soe.tri` | `soe.cgs` | `soe.addA` | `elem_kern` | `mpm.a.prep` | `upd_constr` | `dc.setSz` | `dc.named` | `unbal_o` | `cmt_el` | `cmt_rec` | `residual` | `1t_deep` | `1t_coarse` | `nt_coarse` | `deep_ovh` | `1tc/ntc` | `threaded?` | `max_load` | `rankable?` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| `d1b_brick_concrete3d[smoke]` | 600 | 6 | 0.0% | 1.5% | 1.5% | 0.0% | 0.7% | 64.1% | 0.0% | 31.6% | 0.0% | 0.0% | 0.3% | 0.1% | 0.0% | +0.14% | 0.766s [0.760..0.769] | 0.757s [0.741..0.819] | 0.955s [0.926..0.967] | 1.01x | 0.79x | no (refused -> serial) | 32.0% | NO |
| `d3a_explicit_cdl_cpml[smoke]` | 1944 | 20 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 57.2% | 0.0% | 23.3% | 0.0% | 0.0% | 0.9% | 18.3% | 0.0% | +0.03% | 1.153s [1.139..1.230] | 1.114s [1.099..1.161] | 1.115s [1.075..1.122] | 1.03x | 1.00x | no (refused -> serial) | 51.7% | NO |

**`elem_by_type` breakdown:** `d1b_brick_concrete3d[smoke]`: LadrunoBrick(33002)=64.1% ; `d3a_explicit_cdl_cpml[smoke]`: LadrunoCPML(33023)=55.1%, LadrunoBrick(33002)=2.0%

**Per-deck max host CPU load (threshold <= 10.0%):** `d1b_brick_concrete3d[smoke]`: before=32.0%, after=30.0% ; `d3a_explicit_cdl_cpml[smoke]`: before=38.7%, after=51.7%

**Step-1 one-off costs (excluded from steady-state shares & wall times):** `d1b_brick_concrete3d[smoke]`: 0.0445s ; `d3a_explicit_cdl_cpml[smoke]`: 0.1155s

> **NOT RANKABLE decks (excluded from §4 ranking & §5 picks):** `d1b_brick_concrete3d[smoke]`: cpu_load_before=32.0%>10.0%, cpu_load_after=30.0%>10.0% | `d3a_explicit_cdl_cpml[smoke]`: cpu_load_before=38.7%>10.0%, cpu_load_after=51.7%>10.0%

## 2. Size Trend (`p` at S vs L)

_Only one size present in this run (smoke/partial run); full S vs L table prints once S and L both run._

## 3. Phase -> Package Map & Coverage Gaps

- **`soe.{symbolic,factor,trisolve,cgs}`** -> `260a` (`soe.factor` on Pardiso via opt-in `iparm`; `soe.addA`/`dc.setSize` on MUMPS only), `260g` (mixed-precision IR, unranked), `260m` (AMG, Wave 3).
- **`soe.addA`** -> `260a` (MUMPS only), `260f` (persistent assembly/fast scatter, unranked).
- **`elem_by_type`** -> `260c1` (`LogStrain` finite-strain, N/A in small-strain D1), `260c2` (`LadrunoConcrete3D` on D1b, blocked `#1154`), `260h` (analytic return maps, unranked & blocked `#1153/#1154`), `260k` (Wave 3).
- **Explicit `LadrunoCPML` + `commit.elements` + `ExplicitBathe formTangent` + SMS `soe.trisolve`** -> `260b` (`Unblocked`).
- **`upd_constrained` (`TransformationConstraintHandler::enforceSPs`)** -> unbucketed double element update on constrained elements inside `newStep`/`update`. **Handler dependence (Finding 9):** on D1b this share (`d1b_brick_concrete3d[smoke]`=31.6%) comes from `constraints Transformation` (fixed base + `equalDOF` top cap); a production D1b using a different constraint handler would not see it.
- **`260c2` conservative ranking note (Finding 5):** `_p_260c2` ranks D1b using ONLY the `LadrunoBrick` (`33002`) `elem_by_type` share (conservative), excluding `upd_constrained` (which would raise D1b's total Concrete3D-evaluating share if counted, because `enforceSPs` calls `LadrunoBrick::update()` on constrained elements).
- **`mpm.a.prepare`** -> MPM host work (program `#1075`, tracked cleanly in its own column).
- **`commit.record`** -> recorders (outside WP-260; flagged automatically below if >= 10 %).
- **Coverage gaps >= 10 % (computed from rankable rows):** none.

## 4. Plan §0.5 + Brief §4 Per-Deck Ranking (`p >= 10 %`, `G >= 1.10x`, rule-2 `s`)

| Pkg | Tier | Status / Blocker | Plan Hypothesis | Rule-2 `s` | Source | `d1b_brick_concrete3d[smoke]` | `d3a_explicit_cdl_cpml[smoke]` | Peak L `p / G` | Gate |
|---|---|---|---|---:|---|---: | ---:|---:|:---:|
| `260a` | G-BYTE+G-OPT | Unblocked | 1.3-2x PARDISO factor | 1.15x | plan §1 rule-2: 1+(1.3-1)/2 (MUMPS 3-20x N/A on Pardiso) | N/R | N/R | — | BELOW |
| `260b` | G-BYTE/G-EQUIV | Unblocked | 1.4-2.2x explicit | 1.20x | plan §1 rule-2: 1+(1.4-1)/2 | N/R | N/R | — | BELOW |
| `260c1` | G-BYTE/G-EQUIV | Unblocked | 3-6x finite-strain GP | 2.00x | plan §1 rule-2: 1+(3.0-1)/2 | N/R | N/R | — | BELOW |
| `260c2` | G-BYTE/G-EQUIV | Blocked (#1153/#1154) | 1.5-4x Concrete3D/HS | 1.25x | plan §1 rule-2: 1+(1.5-1)/2 | N/R | N/R | — | BELOW |
| `260d` | G-EQUIV/G-OPT | Unblocked | 2-5x adaptive runs | — | changes step count -> whole-run bench required | N/R | N/R | — | not phase-rankable |
| `260e` | G-BYTE | D4 + build bench | to be measured | — | global allocator/OMP -> whole-run bench required | N/R | N/R | — | not phase-rankable |
| `260f` | G-BYTE/G-EQUIV | Blocked (#1147/#1153) | no number in plan | — | unranked (threading belongs to #1147/#1153/#1154) | N/R | N/R | — | unranked |
| `260g` | G-OPT | Unblocked | no number in plan | — | unranked (mixed-precision PARDISO requires bench) | N/R | N/R | — | unranked |
| `260h` | G-EQUIV/G-OPT | Blocked (#1153/#1154) | no number in plan | — | unranked (SAS-RK / analytic Concrete3D micro-bench needed) | N/R | N/R | — | unranked |
| `260i` | G-OPT | Post-Wave 1 | no number in plan | — | unranked (global IMPL-EX / H5DRM v2) | N/R | N/R | — | unranked |
| `260j` | G-OPT | Post-Wave 1 | no number in plan | — | unranked (contact BVH) | N/R | N/R | — | unranked |
| `260k` | G-OPT | Rule 6 | Wave 3 | — | not eligible for Wave 1 pick | N/R | N/R | — | not eligible (Wave 3) |
| `260l` | G-OPT | Rule 6 | Wave 3 | — | not eligible for Wave 1 pick | N/R | N/R | — | not eligible (Wave 3) |
| `260m` | G-OPT | Rule 6 | Wave 3 | — | not eligible for Wave 1 pick | N/R | N/R | — | not eligible (Wave 3) |
| `260n` | G-OPT | Rule 6 | Wave 3 | — | not eligible for Wave 1 pick | N/R | N/R | — | not eligible (Wave 3) |
| `260o` | G-OPT | Rule 6 | Wave 3 | — | not eligible for Wave 1 pick | N/R | N/R | — | not eligible (Wave 3) |

## 5. Proposed Pick

- No pick: partial table (D9). Missing: `d1a_brick_sanisand[S]`, `d1a_brick_sanisand[L]`, `d1b_brick_concrete3d[S]`, `d1b_brick_concrete3d[L]`, `d2_cpml_wave[S]`, `d2_cpml_wave[L]`, `d3a_explicit_cdl_cpml[S]`, `d3a_explicit_cdl_cpml[L]`, `d3b_explicit_bathe_cpml[S]`, `d3b_explicit_bathe_cpml[L]`, `d3c_explicit_sms_cons[S]`, `d4_mpm_modeA[S]`, `d5_mcft_2cld12[S]`.
