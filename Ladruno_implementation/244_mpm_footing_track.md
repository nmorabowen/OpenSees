# WP-244 — the MPM footing box keeps the soil (fix #1112)

PR #1136, program #1075. Gates pre-registered in `wp-briefs.md` "WP-244" (`de4d8df6c`) before any
footing number on the fix was read. Builds: Windows `ae3e57ad6` for the local runs (the branch changes no
C++; the later merge of `origin/ladruno` brought RC/fiber sources only, none under `ladrunoMPM`), and on
esmeralda `0a8613c01` (the branch, OpenSeesPy) with two BLAS configurations: spack oneMKL 2024
(the WP-238/239 configure) and the default apt MKL 2020 (232j's warrant configuration). Evidence: `tests/data/mpm/wp244/`.

## The defect, made precise

Every footing in the 232f / 232j L3 / L4 / 232k S3 cases is
`gridvel footing box <x0> 0.0 <B/2> 1.0e9 v 0 -0.002 start t0 track y`.

- `track y` lowers the box's lower bound to −s, so the footing prescribes the node rows at or
  above its base; the next row joins only once s reaches it.
- The soil's top particle row (ppc 2) is centred at −s − h/4. As s → h, the gap from the lowest
  prescribed row (y = 0) to the soil grows to a full cell: the soil next to the footing
  interpolates more and more from free nodes, so the footing loses its grip. In the footings
  this shows as a gradual sink (232f with anti-locking off: 143 → 92 kPa over s/h 0.64–0.96) or a
  one-step collapse (HS r3: 868 → 68 kPa in the step 0.752h → 0.768h), and q never reaches zero.
- The extreme case, isolated by P1: once the overlap is thin enough, the prescribed row's nodal
  mass falls below the mode-A mass cutoff and it leaves the active set (`want`,
  `LadrunoMPMImplicit.cpp:754-768`). A free elastic block carried by the box then gets no
  prescription for a step and lags exactly Δs from the step that starts at s = 0.9h; with
  `mpm massCutoff 1e-14` it does not. In the footings the collapse starts before that point
  (the critic's estimate: the rigid-carry cutoff at about 0.84h, the axis node at about 0.78h).
- `translate y` (232k D17) is **not** a fix: these boxes are unbounded above (hi_y = 1e9), so
  `translate y` and `track y` select the same rows.

## The fix

The box's lower bound starts at **−h/2** instead of 0, still `track y`: a node row is prescribed
once the footing base is nearer to it than to the row above. The lowest prescribed row then
always has at least h/2 of soil on its support. No C++ change; `track` and `translate` keep their
meaning. Edited: `tests/mpm_implicit_cases.txt` (232f), `tests/mpm_antilocking_cases.txt` (L3,
L4), `tests/ladruno_mpm_soil_s_harness.py` (`FOOTING_LEAD = 0.5`, S3).

Before s = h/2 the selected rows are unchanged, so every run that ends before s = h/2 is
**bit-identical** on one build at one thread count (checked: S3 HS r1 / r2, the P1 elastic
footing, and every lead-0 / lead-0.5 pair below before h/2). Those records are left as committed.
Thread count alone moves results by ~1e-13 (OpenMP reduction order), so identity holds only at
equal `OMP_NUM_THREADS`.

## Results

"Drop" = the largest one-step fall of q after s = 0.01 m, as a fraction of the running maximum
(the F1 measure). "Drawdown" = the largest fall from the running maximum to any later point; it
also catches a gradual loss.

| Run | Box | Reaches s/h | Worst drop | Drawdown | Failed attempts | Note |
|---|---|---|---|---|---|---|
| P1 free block, h = 0.25 | 0 | 2.5 | — | — | 0 | lags one step (Δs) from s = 0.9h: **strict xfail** |
| P1 free block | −h/2 | 2.5 | — | — | 0 | σ ≤ 2.1e-10, displacement error 2e-15 m: **pass** |
| F1 HS S3 r2, to 2h | 0 | 2.00 | 55.2 % @ 0.97 | 86.5 % | 193 | the #1112 reproduction |
| F1 HS S3 r2, to 2h | −h/2 | 2.00 | 4.0 % @ 1.78 | 5.0 % | 272 | **pass** (gate 5 %) |
| F1 HS S3 r3, to 2h | 0 | 0.91 (stop) | 92.2 % @ 0.77 | 92.2 % | 100 | the #1112 reproduction |
| F1 HS S3 r3, to 2h | −h/2 | 2.00 | 1.5 % @ 1.36 | 1.5 % | 301 | **pass** |
| L4 r2 extended (Windows) | 0 | 0.92 (stop) | 53.0 % @ 0.81 | 80.3 % | 385 | |
| L4 r2 extended (Windows) | −h/2 | 1.01 (stop) | 7.6 % @ 0.26 | 14.5 % | 467 | the 0.26 drop and the drawdown precede h/2: same in both boxes (DP softening) |
| L4 r2 extended (esmeralda) | 0 | 0.87 (stop) | 70.1 % @ 0.81 | 93.0 % | 450 | |
| L4 r2 extended (esmeralda) | −h/2 | 0.99 (stop) | 5.3 % @ 0.78 | 15.0 % after the join (127.8 @ 0.70 → 108.6 @ 0.85) | 495 | post-join and inside the #1112 window; DP softens, so the box and the material cannot be separated here (HS shows ≤ 2.6 % in the same window) |
| 232f footing, anti-locking off, r2 | 0 | 1.57 (stop) | 4.9 % @ 0.98 | 45.8 % | 1759 | q sinks 143 → 92 over s/h 0.64–0.96, then snaps back to 154 at s = h |
| 232f footing, anti-locking off, r2 | −h/2 | 1.29 (stop) | 0.4 % | 4.5 % | 659 | a smooth 156–159 kPa plateau |
| 232f footing, anti-locking off, r3 | both | 0.16 (stop) | — | — | 103 | locking stops it before the box matters (as 232f recorded) |
| S3 SANISAND v0 deck, r2 (F4) | both | 0.03 (stop) | — | — | 90 | material-side stall, box-independent; reported |

**F3, the extension.** Neither L4 r2 run reaches the pre-registered s/h 1.6. The new box stops
at s/h 0.99 (esmeralda) and 1.01 (Windows), against the old box's 0.87 and 0.92. The failed
attempts accumulate evenly, 20–37 per 0.08 s/h. The one burst (+136) is at s/h 0.32–0.40, before
h/2, and it is the same for both boxes. The last step hits the halving floor at s ≈ h, where the
top particles' domains leave row 0's support. HS passes the same point at r2 and r3. The records
hold totals only; the per-interval counts above come from the run logs.

**Residual of the fix.** Two things remain.
- **q steps up when a node row joins** (s/h = k + ½). The soil between the two prescribed rows
  then moves with the footing. The step is +2 to 6 % for HS and 232f, +5 to 8 % for L3 r3, and
  +11 to 16 % for L4 r2 (104.6 → 117.0 kPa at t = 75). For L3 r3 the plateau is 560 kPa,
  against the FE B-bar's 543 (the old box gave 526, then collapsed). L3 r3's reported row
  (q_max within 8 % of (2 + π) k) now reads 1.095; it is a report row, not a gate.
- **q dips when a particle row crosses the lowest prescribed row** (s/h ≈ k + ¾). The dip is
  1.5 to 4 % for HS, which cannot soften. For L4 r2 it falls inside a 15 % post-join drawdown
  that DP softening may share.

The L4 r2 gate value now comes after the join: q_max is 127.77 at t = 98, where it was the
pre-join 120.26. It passes either way.

## Changed records — for the owner's acceptance (bands unchanged)

| Record | Before | After | Gate |
|---|---|---|---|
| `soil/s3_mpm_hs_r3_fbar.json` (S3 HS r3, reported) | q(s/B 0.05) 155.3 kPa = 0.16 × FE | 928.6 kPa = 0.969 × FE n16 (958.7) | reported row; would now sit inside the r1 / r2 ±10 % band |
| `antilocking/footing_j2_r3.json` (L3 r3) | window 500.77 (pass); q_max 525.76, collapse at s/h 0.78 | window 500.77 (pass, unchanged: before the join); q_max 562.84 (reported) | unchanged |
| `antilocking/footing_mildcone_r2_al.json` (L4 r2, esmeralda) | q_max 120.52 (0.922 × FE 130.66), reaches 0.2 m | q_max 127.77 (0.978 × FE), reaches 0.2 m, 402 failed | pass → pass, inside [117.60, 143.73] |
| `antilocking/footing_mildcone_r3_al.json` (L4 r3) | see "#1107" below | | |

## #1107's MPM half

**The configuration matters first.** L4 r3 sits on a stagnating ladder, so where it gives up
depends on round-off from the BLAS build. On esmeralda, the branch's C++ with the default configure
(apt MKL 2020, the sequential preload, as 232j's warrant build) reproduces the committed 232j record
**bit for bit**: every q equal, 502 failed attempts, a stop at t = 57 (s = 0.094 m). With spack
oneMKL 2024 the same case stops at s = 0.010 m, and so does 232j's own merge commit `3e738308c`.
The Windows build stops at 0.010 m too, at either thread count. So nothing regressed, and only
the apt-MKL configuration can answer #1107.

**On that configuration** (`tests/data/mpm/wp244/apt_l4_r3_lead*.json`, build `0a8613c01`):

| Box | Stops at | Failed attempts | q_max (band [110.28, 134.78]) | After the join |
|---|---|---|---|---|
| 0 (the 232j record) | t = 57, s = 0.094 m, **s/h 0.75** | 502 | 115.64 @ s/h 0.14 | — |
| −h/2 | t = 74, s = 0.128 m, **s/h 1.02** | 725 | 115.64 @ s/h 0.14 (unchanged) | 103 → 115 by s/h 0.7; dip to 100.0 at s/h 0.83 (a 13.5 % drawdown from 115.6); 102–103 to the stop |

The two runs are bit-identical before s = h/2. With the fix the run passes s/h 0.75 and stops
at s ≈ h, at a similar failure rate per unit of s (670 per h before the old stop, 820 per h
after it), exactly as L4 r2 does on both platforms.

**Answer.** The 232j stop at s/h 0.75 sat at the onset of the #1112 window, and the fixed box
moves it to s ≈ h. Its cause is still the Newton stagnation of non-associated DP, which #1107
tracks and which the FE shares. L4 r3 stays a miss: it does not reach s = 0.2 m, and q_max is
unchanged. #1107 stands for both halves; only the stop location of its MPM half was set by the
artifact.

**The L4 r3 record.** The committed `antilocking/footing_mildcone_r3_al.json` and the
`L4_R3_MISS` text (stop t = 57) describe the old box. With the fix, on the record's own
configuration, the stop moves to t = 74 (s = 0.128 m). Updating the record and the miss text is a
gate-record change: it is the owner's call. The fixed-box warrant record is prepared on the apt
build and is not committed.
