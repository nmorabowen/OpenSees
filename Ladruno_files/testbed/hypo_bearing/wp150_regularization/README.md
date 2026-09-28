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
| `acoustic_vec.py f1.npz,...` | Plane-strain acoustic tensor of the continuum tangent at every GP, against an associated control and the ADR-90 V4 viscous blend | `out_acoustic.txt` |

The inputs came from the WP-138 analysis copy. The legs are E_B (B/8, SAS-ME, TolR 1e-4) and E_B16 (B/16), Esmeralda
jobs 148583 and 148586 (runs `~/ladruno_wp138/deck/runs/<leg>/`). The checkpoints are not in the repository: they are
3 MB each, and copies live with the WP-138 orchestration. To regenerate them, re-run the leg with `--ckpt-every 5`.

The parameter set is hard-coded in `h_decomp.py`. Pass `NAME=value` after the file list to override one, for example
`A0=0.001` for the ablation S4 set.
