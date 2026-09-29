# WP-150 reproducers: why the WP-138 footing walls, and why it localizes

Memo: `Ladruno_implementation/150_sanisand_regularization_memo.md`.

The scripts are post-processors. They read the WP-138 Esmeralda checkpoints (`field_*.npz`, written by
`footing_ab.py`: committed stress, the 26-slot `state`, ψ, the SAS counters and GP coordinates) and mirror the
SANISAND kernel formulas (`LadrunoSANISANDSasME.cpp:335-411`, `ManzariDafalias::GetStateDependent` /
`GetElasticModuli`, e_G = e_init, p_r = 0, the campaign parameter set). No engine is loaded.

| script | question | output |
|---|---|---|
| `h_decomp.py f1.npz,f2.npz,...` | Is any committed GP near the SAS-ME refusal H = Kp + 2G − K·D·qv ≤ 0? It prints the term split and the lowest-H GPs. | `out_h_decomp.txt` |
| `refuser_stats.py <analysis dir>` | H split and cumulative SAS counters (re-seats, reversal rejections) at the GPs that refused at the wall (`tables/floor_refusers_in_band.csv`) | `out_refusers.txt` |
| `t5_element_physics.py [eps_max] [SET=toyoura] [E0=..] [NAME=value]` | T5: drained plane-strain and triaxial compression on the EXACT WP-134 oracle (no engine), with Bolton (1986) stress–dilatancy checks | `out_t5*.md`, `out_t5*.json` |
| `t5_pdmy_control.py <dist/bin>` | T5 control: the fork's WP-133 PDMY03 stand-in (NOT TIMs' PDMY01), drained plane strain | `out_t5_pdmy03_standin.md` |
| `t6_capacity_bands.py` | T6: the classical rough-strip capacity band of the deck (Martin 2005 exact N_γ, exact N_q) | `out_t6_capacity.md` |
| `lode_split.py f1.npz,...` | The §2.1 non-elliptic GPs split by the Lode angle of n (does c < 7/9 extension non-convexity drive the bands?) | `out_lode_split.txt` |
| `r2_analysis.py LABEL=RUNDIR ... --family b4,b8,b16 --h b4=0.375,... [--png]` | R2 (memo §9.1 steps 3–6): q at matched s/B, mesh-family differences, contraction and a Richardson limit; band path, inclination, FWHM/h and w2 under the right footing edge | stdout (markdown), optional PNG |
| `acoustic_vec.py f1.npz,...` | Plane-strain acoustic tensor of the continuum tangent at every GP, against an associated control and the ADR-90 V4 viscous blend | `out_acoustic.txt` |

The inputs came from the WP-138 analysis copy. The legs are E_B (B/8, SAS-ME, TolR 1e-4) and E_B16 (B/16), Esmeralda
jobs 148583 and 148586 (runs `~/ladruno_wp138/deck/runs/<leg>/`). The checkpoints are not in the repository: they are
3 MB each, and copies live with the WP-138 orchestration. To regenerate them, re-run the leg with `--ckpt-every 5`.

The parameter set is hard-coded in `h_decomp.py`. Pass `NAME=value` after the file list to override one, for example
`A0=0.001` for the ablation S4 set.

## R2 deck support: B/4 and mesh-orientation variants

`footing_ab_meshperturb.patch` applies to the WP-138 Esmeralda deck `footing_ab.py` (base sha1
`6425ff3ca7c38b6a8808ca80979b41e925448c8a`, the orchestrator's copy). Put `mesh_perturb.py` next to it.
- `--mesh b4`: B/4 fine band, 46 × 14 = 644 elements, 5 footprint nodes.
  - Graded counts 7 (x) and 8 (y), with 6 fine rows.
  - The node-count and footprint assertions are generalised; b8/b16 are unchanged.
- `--mesh-perturb shear:DEG | jitter:A[:SEED[:KEEP]]`: moves interior fine-zone nodes only (see the `mesh_perturb.py`
  docstring). The Gauss-point coordinates and areas logged and saved use the true distorted geometry.

**Validated locally** (DP 38°, engine dd107e5aa copy, 3 push steps each):

| variant | worst min/max detJ | 1-D K0 patch max rel err after gravity |
|---|---|---|
| none (b8) | 1 | 2.1e-12 |
| **shear:15** (the primary orientation leg) | 0.974 | 8.8e-3 |
| jitter:0.1 (KEEP 1; the secondary leg) | 0.546 | 4.6e-2 |
| jitter:0.2 (KEEP 1) | 0.239 | 9.0e-2 |
| jitter:0.2 KEEP 0 (avoid) | 0.239 | 2.6e-1 |
| b4 | 1 | 5.3e-13 |
| b4 shear:15 | 0.950 | 1.7e-2 |

- Base reaction is exact in every case.
- The K0 patch error is a discretization error: distorted bilinear quads do not carry the 1-D self-weight field
  element by element. It is largest where the geostatic stress is smallest (near the top of the fine band).
- The **default path is byte-identical** to the unpatched deck: b8, 3 steps, steps.csv field-for-field apart from wall
  time.

## GATE 0: the DM04 Toyoura reference sand (memo §12)

| script | what | output |
|---|---|---|
| `gate0_toyoura_oracle.py` | DM04 Figs. 5–9 test matrix (17 triaxial tests) on the exact oracle, `paper` and `uw_model` options | `out_gate0_oracle.json` (not committed; ~20 s) + stdout |
| `gate0_cxx_driver.py <bin> <site> <spec> <out>` | the C++ LadrunoSANISAND (SAS-ME) through a mixed-control `ladrunoSANISANDReplay` loop (TXu / TXd / PSd), variants e.g. R1 off/on; run with `python -S` | `out_gate0_cxx.json`, `out_t5_cxx_*.json` (not committed) |
| `gate0_compare.py <oracle.json> <cxx.json>` | max \|Δq\|/q_max, C++ vs oracle and R1 on vs off | `out_gate0_compare.md` |
| `gate0_overlay.py <dafalias2004.pdf> <json> <sets> <outdir>` | overlays on DM04's own figure panels, from the reader's copy (gridline-calibrated axes). The images are NOT committed (copyright) | PNGs in `<outdir>` |

The specs are `gate0_cxx_spec.json` and `t5_cxx_spec_toyoura_e0.643.json`. The large JSON outputs are regenerable and
kept out of git.
