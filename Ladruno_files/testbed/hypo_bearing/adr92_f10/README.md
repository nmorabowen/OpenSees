# adr92_f10 — the IMPL-EX self-weight wall (WP-102 / ADR-92 F10)

The test bed behind
[`Ladruno_implementation/92b_implex_selfweight_wall_note.md`](../../../../Ladruno_implementation/92b_implex_selfweight_wall_note.md).

**A diagnosis, not a feature.** Nothing in `SRC/` was touched or is needed.

## The question

TIMs measured, on a plane-strain strip footing on SELF-WEIGHT SANISAND with
`-implex -implexControl`, that the control's refusal counter climbs from the
first push step (204 120 at step 3 on 9 720 Gauss points), that every refusal
cuts the step, and that the harness step falls below its floor at `s/B = 0.0125`
— while the same material class and the same build carried the weightless ADR-95
slab campaign to `s/B = 0.15`. Minimum `p'` in the block is 3.8 kPa, so this is
not the ADR-93 apex/zero-confinement wall.

Which `-implexControl` criterion fires, at which Gauss points, in what state —
and is there a flag or a staging that carries a self-weight strip to a peak?

**Answer (leg N):** remove `-implexControl`. Bare `-implex` with the same
doubling controller reaches the target `s/B = 0.05` in 104 steps, 0
subdivisions, 58 s, refusal ledger `0/0/0/0` — which is how the ADR-95 campaign
that reached `s/B = 0.15` was run. With the control on, the refusal COUNT is set
by the harness's growth rule and the step FLOOR is set by the control's own
`implexPrimed` bare `> 0.0` test. Full verdict and tables in the note.

## The deck

| | |
|---|---|
| geometry | `B = 1.5` m strip, box `15B x 12B` (22.5 x 18 m), one element thick, `u_y = 0` everywhere (plane strain) |
| mesh | graded, `h0 = 0.5` m under and near the footing, `r = 1.35` outward — **285 `LadrunoBrick -formulation bbar`, 640 nodes, 1 920 DOF, 2 280 Gauss points** |
| soil | `LadrunoSANISAND`, the ADR-92 CP1 / ADR-86 §5 Gorini set (`tests/test_ladruno_sanisand.py::_PARAMS`, `e_init = 0.6944`), `IntScheme 1 TanType 2 JacoType 1 TolF = TolR = 1e-7`, `-Presidual 0 -Pmin 0.0101 -maxSubsteps 1000` |
| loads | `gamma' = 9.81` kN/m3 through `eleLoad -type -selfWeight` against the element's `-b`, 7.65 kPa surcharge OUTSIDE the footprint, 18.4 kN/m of footing dead load spread over the footprint |
| K0 | the **nu\* device**: `nu* = K0/(1+K0)` passed as the material's Poisson ratio, so the stage-0 elastic gravity state is exactly `K0` at every depth |
| footing | rigid and ROUGH (`u_x = 0` on the footing nodes), the ADR-92 CP1 convention |
| staging | confine (stage 0) -> `updateMaterialStage 1` -> push |
| push | `LoadControl(-ds)` on an `sp` pattern under a `Linear` series with the `Transformation` handler — the idiom the guide's warning box makes a precondition for `-implexControl`; **never `DisplacementControl`** |
| controller | R3's adaptive halve/double with `DS_BASE = 2e-5`, `DS_MIN = 2e-7`, `DS_MAX = 1e-3`, `GROW_AFTER = 6`, `SUBDIV_BUDGET = 80` (pinned, see `tests/test_r3_prandtl_collapse_gate.py` CONSTRAINT 2) |
| controls | gravity resultant identity, the 1-D geostatic patch (`3.1e-13` measured), `eta_max/M_c < 1` at the flip (`0.7275` measured) |

The architecture is `sanisand_tau0_band.py`'s (ADR-90 WP-A2), which is R3's;
`_graded` / `_n_graded` are transcribed from R3 rather than imported, because
that module imports `_testbed` and would bind a different engine.

**Declared differences from the TIMs deck** (state them beside any comparison):

* the fork's `LadrunoSANISAND` takes `nu` as a positional constant with no way
  to put it back after the K0 stage, so `nu*` is **held for the whole leg**, not
  temporary. At `K0 = 0.455` this is nearly free — `nu* = 0.31271` against the
  material's own calibrated `nu = 0.3129` — so leg B is simultaneously the
  "K0 reached natively" control. At `K0 = 0.818` (`nu* = 0.45`) it is not;
* the calibration is the fork's CP1 set, not TIMs'. Measured `psi` on this deck
  is `-0.128 … -0.109` against TIMs' quoted `~ -0.13`;
* 2 280 Gauss points against TIMs' 9 720, and `B = 1.5` m on a coarser mesh;
* the footing is driven by a prescribed settlement on its nodes, not through a
  `LadrunoKinematicCoupling` master node.

## Engine

The **pinned release build of `ladruno` tip `9c2f964`** (`ladrunoBuild()` =
`9c2f964eae3bbd1a055c3ede81381a6c601a982b`). Nothing was built for this WP.
Override with `LADRUNO_DIST_BIN`; `LADRUNO_F10_EXPECT_BUILD` pins the hash and
accepts `any` (and then says so loudly).

## Files

* `f10_selfweight_wall.py` — the driver. Three commands:
  * `leg` — a full adaptive push, one row of the leg table;
  * `census` — a few push steps at a FIXED `ds` with the control tolerance set
    so high that nothing can refuse, then every Gauss point's `implexDetail`
    read back beside its `p'`, `eta`, `psi`, `M^d`, depth and distance from the
    footing edge. This is the measurement a walled leg cannot make;
  * `probe` — walk to a fixed settlement on a refusal-free path, then take ONE
    step of a given size and census it. Run once per `ds` from separate
    processes, this is a controlled **step-size refinement at a fixed state**:
    does the bulk error field scale with the step, or not? (It does — first
    order.) **Two limits, stated because the probe was over-read in the first
    cut:** it measures the BULK field and reports zero over-tolerance points at
    every `ds`, so it cannot see the handful of seizing points the leg's own
    throttled warning lines carry; and its maximum-error point runs at `f = 0`
    at every `ds` (guarded elastic-predictor drift, not extrapolation) while `f`
    itself scales with `ds` through the sweep, so `O(ds)` and `O(f)` are not
    separated by it.
* `f10_summary.py` — reduces the JSONs and CSVs to the note's tables. Imports no
  engine, so an old campaign's artefacts can be re-reduced on any box.
* `out/` — the measured artefacts (per-leg CSV + JSON, per-census CSV + JSON,
  per-probe JSON), the run scripts, and `out/refusal_warnings_<leg>.txt` — the
  throttled `-implexControl REFUSES` lines, extracted and committed because they
  are the evidence for the note's §4 (the `dt`-INDEPENDENT refusal family that
  the error-field probe cannot see). **Raw engine logs are not committed**
  (`.gitignore` here): tens of MB of throttled warning text plus a banner per
  process, reproduced exactly by re-running the scripts. `out/summary.md` carries
  the reduced tables that the note quotes.

## Running

```sh
cd Ladruno_files/testbed/hypo_bearing/adr92_f10
python3.12 -u f10_selfweight_wall.py leg    --leg B --h0 0.5 --wall 400 --out out
python3.12 -u f10_selfweight_wall.py census --leg B --h0 0.5 --census-ds 4e-5 --out out
python3.12 -u f10_selfweight_wall.py probe  --leg B --h0 0.5 --probe-s 0.0085 --census-ds 1e-5 --out out
python3.12    f10_summary.py out
```

`out/run_legs*.sh`, `out/run_census.sh` and `out/run_probe.sh` are the campaign's
own batch scripts, committed so the reported set is reproducible verbatim.
Measured wall time for the whole campaign: **about 84 minutes** across 17 legs,
16 censuses and 7 probes on the reference box (Windows 11, Intel oneAPI build,
`system Pardiso -matrixType 0`). Four legs (D, G, K, L) terminated on their
wall-clock budget rather than on the physics, and **the box was not attested
idle** — treat every wall number as an upper bound, and do not compare wall
times across legs that ran at different times.

**One writer per leg.** The driver refuses to open `out/f10_<leg>.csv` if another
process touched it in the last 180 s (`F10_FORCE=1` overrides) — the same guard
the ADR-79 runner carries, added here after an accidental double launch
interleaved rows into legs L and M and tore a line. The physics was unaffected
(leg M reproduced `s/B = 0.011286` to the digit on the single-process re-run) but
both wall times were wrong. `f10_summary.py` also drops torn lines.

## The legs

| leg | what it changes against B |
|---|---|
| **N** | **`-implexControl` REMOVED** — bare `-implex`, growth factor still 2.0. THE DECISIVE ARM, and the way the ADR-95 reference campaign ran it |
| **N1** | N with the growth factor also pinned at 1.0 — separates "no control" from "no growth" |
| **A** | weightless, uniform 10 kPa surcharge — the ADR-95 campaign's condition |
| **B** | **the reported configuration**: self-weight, `K0 = 0.455` via `nu*`, `-implexControl 0.05 0.01` |
| **C** | `K0 = 0.818` via `nu* = 0.45` — the ADR-95 slab's own Poisson ratio; CONTRACTANT at rest |
| **D** | `-implexFactor controlIter` |
| **E** | a hold at the flip + a 2e-6 m first push step, growth factor 1.5 |
| **F1** | `tol = 0.1` (the C++ default since WP-92d) |
| **F2** | `reductionLimit = 0.5` — the reduction floor raised so it can actually bind |
| **F3** | `tol = 0.5` |
| **G** | implicit (no `-implex`) — the reference reach |
| **H** | a heavy 100 kPa uniform surcharge on top of self weight (min `p'` up two decades) |
| **I** | `-implexGuard off` |
| **J** | `-implexGuard off` **and** `-implexTrialGuard off` |
| **K** | the controller's growth factor pinned at **1.0** — `ds` never doubles, so the clock ratio `f` can never exceed 1 |

The verdict, the tables and the recommendation are in the note.
