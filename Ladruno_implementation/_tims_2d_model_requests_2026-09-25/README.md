# Attachments to `_tims_2d_model_requests_2026-09-25.md`

`ring_points_b8.csv`, `ring_points_b16.csv` — material states of `LadrunoSANISAND` Gauss points in the low-confinement ring beside a rigid strip footing, for replay at a single material point (F18, F21). Forty rows each: the points with p' < 10 kPa, sorted by η/M^b descending.

Source: the TIMs `2d-model` act's explicit central-difference lane (`CentralDifferenceLadruno`), Esmeralda jobs 147068 (B/8 quad, 9 720 points, failed at s/B 0.0207) and 147069 (B/16 quad, 28 896 points, failed at s/B 0.0127), build `48c0e99bc8e28bbb4fdf965015f285f01e16e90d`, `-maxSubsteps 20000`. State read with the model still in memory right after `analyze` returned −2. Material: the campaign set (G0 264.32, nu 0.312885 [the Jaky K0 substitution], e_init 0.6944, Mc 1.3309, c 0.71, lambda_c 0.027, e0 0.83, ksi 0.45, Patm 101, m 0.005, h0 1.3, ch 0.968, nb 3.5, A0 0.05, nd 5.75, zmax 12.5, cz 1100, Den 2.0), `IntScheme 1`, `TanType 0`, `-flipAlphaIn init`, `-Pmin 0.0101`, `-Presidual 0`.

Columns:

| column | meaning |
|---|---|
| element, gp | element tag and Gauss point (1-4) in the act's mesh |
| x_m, y_m | element centroid (m); the footing spans x in [-0.75, 0.75] at y = 0, soil below |
| p_kPa | p' = -tr(sigma)/3 |
| eta | q/p' |
| eta_over_Mb_compression | eta / (Mc exp(-nb psi)), with g(theta) = 1: a screening ratio, not the model's own M^b at the point's Lode angle |
| e, psi | void ratio and state parameter psi = e - e_c(p') |
| sigma_0..5 | stress, Voigt order xx, yy, zz, xy, yz, zx, kPa, OpenSees sign (compression negative) |
| alpha_0..5 | back-stress ratio alpha, same order |
| alpha_in_0..5 | alpha at the last load reversal |
| z_0..5 | fabric tensor |

The `substeps` response read 0 at every point in this dump (it reports the last call only and the failed step had been reverted); that is part of F20(a).
