# `implex_ON` — pre-fix arm (superseded)

This leg (`h1.0_e0.6944`, build `80e5...` / `80e4f5630`) was run **before**
`3c788778f` fixed D2's negative-pseudo-dt gate. Every subdivision rung of
every step refused at the DS_MIN floor (`s/B = 0.00000`, 0 converged steps),
so the push never started — see `run.log` and `a2_h1.0_e0.6944.FAILED.txt`.

The removed `a2_h1.0_e0.6944_engine.log` was 12.4 MB / 34,396 lines, of
which 34,272 lines were one repeated D2 refusal warning (4,896 steps x 7
subdivision levels). It carried no evidentiary value beyond the count above
and was deleted for size (red/blue review F12.3). The post-fix arms'
engine logs (`implex_ctl/`, `implex_noctl/`) are 394 and 8 lines.
