# ADR-95 cross-check #3 — does LadrunoSANISAND wall at the same trigger?

**PREDICTION (written before running).** Note 95 §3b's trigger — first footing-edge GPs
reaching mean stress >= 0 (p -> 0, SANISAND's own compression-positive convention) — is
implementation-independent; SANISAND has no cohesion and its moduli/state (`M^b`, `M^d`,
elastic `G`) go as `sqrt(p)`, undefined at `p <= 0`. Expected: `h20uri` walls at essentially
the SAME station as UW/ASD, as a material REFUSAL/substep exhaustion, not a tangent event;
the linear `h8bbar` control reaches further. **Confirmed only in direction, not mechanism —
see VERDICT.**

## VERDICT

| leg | mode | s/B end | q_end (kPa) | I1_max @ end | reached p>=0? | nsub/nfail |
|---|---|---|---|---|---|---|
| h8bbar (linear, implicit, cap 1000) | **WALL** | 0.0317 | 245.4 | −11.0 | **no** | 32/564 |
| h20uri (quad, implicit, cap 1000) | **WALL** | 0.0130 | 114.8 | −11.7 | **no** | 88/1570 |
| h20uri (quad, IMPL-EX, cap 1000) | **TARGET** | 0.1500 | 2163.9 | **−0.303 (clamped)** | no (held at `-Pmin`) | 0/0 |

**Neither implicit leg reaches the literal trigger.** Both die of WALL-CLOCK (not
FLOOR/BUDGET/refusal-cascade), `n_tension=0` throughout, `I1_max` plateauing ~−11 kPa, an
order of magnitude short of zero. `h20uri` is the more fragile leg (12x less depth, ~2.8x
more failed rungs, 194 vs 142 `Domain::update` failures in the same 5400 s) — note 95's
linear-vs-quadratic split, but via cost, not a clean event. **IMPL-EX is the only leg that
reaches p->0**, and there the ADR-86 `-Pmin` clamp — not a refusal/NaN/corner branch — holds
the footing-edge GPs at p=+0.101 kPa (`I1=3x(-Pmin)=-0.303`) for the rest of the push;
`n_tension` stays 0 by construction. The signature is a saturating `implexRefusals`/
`companion` counter (0->1620 over s/B 0.031-0.150), invisible without that response.

**Methodology deviation, load-bearing.** Deck default `-maxSubsteps 0` (uncapped) made a
single `analyze()` near s/B~0.012-0.025 take 15-45+ min with the CPU actively computing
(`Get-Process` CPU-time deltas confirmed it, not deadlocked) — a blocking call `tmax` cannot
interrupt. Switched to `-maxSubsteps 1000`, the ADR-92 CP1 campaign's own cap, before all
three legs.

## Deck

Geometry/mesh/surcharge/push reused **by import** from `h20_prandtl.py`: graded strip,
`B_FOOT=2 m`, `XLIM=30`, `ZBOT=-20`, `Q0=10 kPa` consistent Q8/Q4 surcharge (whole top,
footing included, reaction-corrected), weightless (`GAMMA=0`), ADR-63 D16 ladders. Footing
**ROUGH** (u_x=0 on footing nodes, ADR-92 convention; `h20_prandtl`'s own footing is smooth).
Material `LadrunoSANISAND`, Gorini's calibrated `_PARAMS` (`tests/test_ladruno_sanisand.py`,
ADR-86 sec.5, `e_init=0.6944`); optionals `IntScheme=1 TanType=2 JacoType=1 TolF=TolR=1e-7`
(ADR-92 deck's choice); flags `-Presidual 0.0 -Pmin 1e-3*P_atm -honorTolR 0 -maxSubsteps 1000`,
`-implex` for leg 3. Staging: `-stage 0` (elastic) under the surcharge ramp, `eta_max/M_c<1`
asserted at the flip, `-stage 1` (plastic) before the push. `q_exact`: **none quoted** —
SANISAND has never produced a capacity in ADR-92 CP1 (`_adr92_cp1_surcharge_results.md`:
every leg WALL/FLOOR, no plateau/peak); this deck's domain and confinement (surcharge, not
self weight) differ from CP1 too, so no ratio is comparable. Census (`count_tension_gps`,
copied from `asd_path_diag.py`) plus per-GP `substeps`->`[substepsTakenInME, capHit]` and,
IMPL-EX only, process-wide `implexRefusals`.

## Exact CLIs

```bash
cd Ladruno_files/testbed/hypo_bearing
export ADR95_DIST=<staged bin, ladrunoBuild()=c945f9a8b5e92783109efcfd25813e6591bbc236>
py -3.12 sanisand_path_diag.py --elem h8bbar --sfrac 0.15 --tmax 5400 --maxsubsteps 1000 --suffix _sanisand
py -3.12 sanisand_path_diag.py --elem h20uri --sfrac 0.15 --tmax 5400 --maxsubsteps 1000 --suffix _sanisand
py -3.12 sanisand_path_diag.py --elem h20uri --sfrac 0.15 --tmax 5400 --maxsubsteps 1000 --implex --suffix _sanisand_implex
```

Outputs: `sanisand_{tag}.csv` / `_tension.csv` per leg, engine logs `sanisand_leg{1,2,3}_*.log`.
