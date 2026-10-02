# RED TEAM — Kimura 1985 Fig. 9 V85.6 vs W_TYR_fig9_856 (2026-09-29)

Scratch files here: k_p*.txt / d_p*.txt (PDF text), fig9_rt.png (own 500-dpi render, not for sharing), rt_dark.npy, rt_calib.npy, p5/p6/p8/p9.png (page renders, not for sharing).

## Table

| # | Item | Status | Evidence | Impact |
|---|---|---|---|---|
| A1 | q axis "x100 kN/m2" -> kPa x100 | OK | Fig. 9 axis label (page 8 render); Table 3 (p.42) observed N_g 174 for V75.1 <-> digitized peak 1236 kPa = 174*0.5*15.9*0.9 = 1245 (0.7 %) | none |
| A2 | S (model mm)/30 -> s/B; B_proto 0.9 m | OK | caption "B = 30 mm, 30g"; CSV header; s_proto = S*30 | none |
| A3 | Stress similarity: 1g prototype, gamma 15.9, B 0.9 | OK | footing_ab.py:616-636; log line 7 | see C7 (g gradient) |
| A4 | Model q: gross reaction/B, full B, full-width model, per m | OK | footing_ab.py:784-786 (W_FOOT = 0), q0 = 0 (log:41); equilibrium check: integral of -s_yy over the top element row and rows at z 0.5/0.96/1.29 m equals q*B + gamma*z*W to 1e-4 at steps 100/1000/2205/2415 | none |
| A5 | Test q includes footing self-weight at 30g; model has none | OK (<1 %) | not stated in paper; est. Al footing 30x100x~20 mm at 30g ~ 16 kPa | <1 % on peak; may matter for the s/B 0.02 secant (see E5) |
| A6 | s measured at footing (REF node coupled to 9 footprint nodes), from post-gravity position; no seating | OK | footing_ab.py:775-779, 811, 855-856; uy0 = -3.96e-3 m (log:40) | none |
| B1 | Axis calibration | OK | own 500-dpi render: x ticks 176/441/717/992/1263.5, y ticks 1507/1228.5/955.5/680/398.5/131.5 -> spacings agree with CSV header to <=1 px equiv. Shear: I fit 0.06 deg vs their 0.4 deg (their fit likely cleaner: mine includes label rows); worst-case S bias 0.06 mm (s/B 0.002) at q 20 | <= 0.002 s/B, <= 30 kPa |
| B2 | Curve identity | OK | label leaders on page render: "V 85.6 %" attaches to the curve peaking 19.5 at S~2.7; "H 86.7 %" to the one peaking 17.3 at S~4.9 | none |
| B3 | Independent re-read, V85.6 (16 pts) | OK | S=0.6:7.50 (CSV 7.44); 1.0:10.79 (10.85); 1.5:14.44 (14.65); 2.0:17.28 (17.51); 2.5:19.19 (19.32); 2.72:19.52 (19.53); 3.0:18.38 (18.06, leader crossing); 3.5:16.01 (15.87); 4.0:14.54 (14.42); 4.5:13.42 (13.34); 4.9:12.63 (12.59); 5.5:11.63 (11.59); 6.0:10.76 (10.72); 6.5:10.00 (9.99); 7.0:9.60 (9.58). Peak 1952 kPa at S 2.72-2.77 mm (s/B 0.091-0.092) vs CSV 1952.7 at 0.0906 | max 30 kPa (1.5 %) at leader crossings, typically 10 kPa |
| B4 | Independent re-read, H86.7 (12 pts) | OK | 2.0:8.27 (8.36); 2.5:10.40 (10.54); 3.0:12.38 (12.59); 3.5:14.38 (14.49); 4.0:15.79 (15.93); 4.5:16.88 (16.95); 4.9:17.33 (17.33, peak); 5.5:16.02 (15.76, leader); 6.0:14.29 (14.25); 6.5:13.17 (13.08); 7.0:12.31 (12.22); 7.5:11.61 (11.53) | <= 25 kPa |
| B5 | Claim numbers | OK | model peak 2015.35 at s/B 0.16721 (steps.csv row 2206); q(0.02) = 328.0 vs test 743.8 -> 1-0.441 = 55.9 % | none |
| C1 | Full vs half model | OK | full width, no symmetry (footing_ab.py:15-16, 266-270) | none |
| C2 | Footing width in mesh = 0.9 m | OK | B_FOOT from --bproto (:616); foot = nodes with |x| <= 0.45, 9 nodes asserted (:266-268); log:23 "fine h = B/8 = 0.11250 m"; summary toyoura.B 0.9 | none |
| C3 | Rigid, guided (u_x, rot fixed), fully rough (both DOFs coupled) | OK model / UNVERIFIABLE test | footing_ab.py:264-270. Paper states rough (glued sand) only for the roughness series (p.37) and the de Beer series (p.36); the anisotropy series' footing roughness is NOT stated (would be in Kimura et al. 1979). Driver docstring asserts "ROUGH" | if the test footing were smooth: rough/smooth = 1.25 at B_N 1.2 m, >1.3 at 0.6 m (Fig. 6) -> model would be 25-35 % high |
| C4 | Embedment 0 | UNVERIFIABLE (likely OK) | not stated for Fig. 9; embedment is a separate section (p.42) | none if surface |
| C5 | Domain 16.67B x 10B = container 0.5 x 0.3 m; sand depth not stated | UNVERIFIABLE | Table 1 p.34; footing_ab.py:568-571 assumes full depth; sides rollers, base fixed | < 2 % on q_peak if sand depth >= 6B |
| C6 | Plane strain vs 100 mm-wide box (L/B = 3.3, steel walls, no lubrication mentioned) | UNVERIFIABLE | Table 1 (width 0.1 m); p.35 "sand 100 mm thick in a 9 mm steel-plated box" | side-wall friction typically raises measured q by 5-15 % |
| C7 | 30g uniform vs radial field (R_eff 1.18 m, depth 0.3 m: +-13 % top to bottom) | UNVERIFIABLE | Table 1; reference radius of "30g" not stated | +-4-8 % equivalent gamma in the failure zone depending on reference depth |
| C8 | Displacement control both | OK | geared-motor loading (Fig. 1); model sp push | none |
| D1 | Argument order | OK | LadrunoSANISAND.cpp:63-66: G0 nu e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz Rho <IntScheme TanType JacoType TolF TolR> flags; numData 18 (:294) | none |
| D2 | nu = 0.3333 in matdesc (LEAD) | OK, deliberate | nu* = K0/(1+K0) gives K0 = 0.5 in the elastic gravity stage (footing_ab.py:628-629, 579-582); after updateMaterialStage 1, setParameter poissonRatio -> 0.05 (:767-771; ManzariDafalias.cpp:945, 985-987 sets m_nu); verified from the elastic tangent C12/C11 = 0.052632 -> nu 0.050000 at push step 1 (log.log:39, 42-43). Push runs DM04 nu 0.05 | none. (Had it stayed 0.333, K = 2.67G instead of 0.78G: stiffer volumetric response) |
| D3 | G0 125, Mc 1.25, c 0.712, lc 0.019, e0 0.934, xi 0.7, m 0.01, h0 7.05, ch 0.968, nb 1.1, A0 0.704, nd 3.5, zmax 4, cz 600 | OK | DM04 Table 1 (d_p5.txt:124-163) | none |
| D4 | P_atm 100 vs 101.3 | OK (choice) | not tabulated in DM04; G ~ sqrt(p*P_atm) -> 0.65 % | < 1 % |
| D5 | Rho 1.6208 t/m3 = 15.9/9.81 | OK | mass only (static run) | none |
| D6 | e_init 0.65172 = 0.977 - 0.856*(0.977-0.597) | OK arithmetic / UNVERIFIABLE source | e_max/e_min are Verdugo & Ishihara's batch; Kimura 1985 gives no e_max/e_min/G_s (grep of all pages). Other Toyoura batches: e_min 0.597-0.616 | with 0.977/0.605: e_init 0.6586 (+0.007) -> d ln(alpha_b) = -1.1*0.007 = -0.8 % -> dphi ~ -0.35 deg -> N_g -7 to -8 % (memo s.14: +22 %/deg at 44-45 deg). +-5-8 % on q_peak |
| D7 | gamma 15.9 with e 0.6517 | ERROR (minor) | 15.9 is the vibro-compacted de Beer series (p.36), not the pouring series; Gs 2.65: gamma_d(e 0.6517) = 15.74; 15.9 <-> e 0.635 (D_r 0.90) | q ~ gamma: model ~ +1 % |
| D8 | Dry sand both | OK | pouring, radiographs; model has no pore pressure | none |
| D9 | K0 = 0.5 (driver default) | UNVERIFIABLE (choice) | footing_ab.py:628; Jaky at phi_cs 31.2 deg. Air-pluviated dense sand: 0.35-0.45 | lower K0 -> lower p', higher initial eta (0.96 at 0.4) -> model softer still; few % on peak, ~10 % on early secant |
| D10 | Initial alpha = r(K0), alpha_in := alpha, z = 0; K0 patch 1.5e-12; sigma_v = gamma z centre/far | OK | log:16, 24-38; Elastic2Plastic at flip | none |
| D11 | Pmin 0.0101, p_r 0, tension cutoff 0.5/1.0 kPa | OK | sibling B 1.2 legs: cutoff 0.5/1.0 -> q_max 2508.8; 0.25/0.5 -> 2513.1 (0.2 %) | 0.2 % |
| D12 | R1 devices (hFloor 1, reseatHyst 1, softCap 0.5) = DM04 variant | UNVERIFIABLE | summary sas_totals: 4.2e8 h-floored, 4.2e6 soft-capped substeps, 9.9e6 reseats; R1-OFF fig9 leg floored at s/B 0.020 (ladruno_r1 TYR_fig9_856_off), so no OFF reference exists at B 0.9 | unknown on peak/post-peak; Gate 0 element-level <= 6e-4 only |
| D13 | Mesh dependence of q_peak | ERROR (for the "+3 %" wording) | B 1.2 m siblings: b8 2508.8 vs b8 shear:15 2729.1 (+8.8 %); b16 still running (s/B 0.06, 1127 kPa) | model's own mesh uncertainty >= 9 % > 3 % |
| D14 | DM04 calibrated D_r 18-64 %, p 100-3000; footing at D_r 86 %, p' 1-100 | UNVERIFIABLE (extrapolation, acknowledged) | memo s.12 | unknown |
| E1 | "Typical" curves; test scatter | OK/note | Fig. 4 dense band at B_N 0.9 m: N_g ~ 280-330 (+-8 %); Fig. 12 lines are fits | +3 % is inside test scatter |
| E2 | s/B at peak in Kimura's own dense tests | note | 0.091 (Fig. 9 V), 0.163 (Fig. 9 H), 0.16 (p.37: S 4.8 mm at peak, B_N 1.2 m), ~0.09 (Fig. 6b rough) | the 0.167 vs 0.091 gap equals the spread within the paper |
| E3 | Post-peak "too ductile" | note | at equal s/B 0.2: 1748 vs 1072; at delta s/B 0.033 past each peak: model 87 % of peak, test 82 % | mostly the late peak, not the softening rate |
| E4 | V vs H | note | model isotropic; DM04 calibrated on specimens loaded normal to bedding = V. "Between V and H" is not a merit | none |
| E5 | s/B 0.02 secant contamination | OK, but note | q0 = 0; secants 12.6 (0.001), 15.1 (0.01), 16.4 (0.02), 17.8 (0.05) MPa per unit s/B: model is concave-UP at the start (stiffens), test concave-down. Surface GPs at p' 0.6 kPa (G ~ 3 MPa); no surcharge, no footing weight | the 56 % is robust in sign; the magnitude depends on the free-surface treatment (footing weight/seating ~ 16 kPa, K0) |

## Ranked errors that change conclusions

1. "Reproduces peak within +3 %" is over-stated. Independent uncertainties each exceed 3 %: model mesh orientation >= 9 % (D13); e_max/e_min choice +-5-8 % (D6); test scatter +-8 % (E1); side-wall friction and g-gradient of the 100 mm box, unquantified 5-15 % (C6, C7); gamma +1 % (D7). Correct wording: consistent with the test within the combined uncertainty (~ +-15 %).
2. Footing roughness of the anisotropy series is unverified (C3). If smooth, the model is 25-35 % high and the peak agreement is spurious.
3. gamma 15.9 belongs to the vibro-compacted series and is 1 % inconsistent with e_init (D7): fix to 15.7 or derive e from gamma.
4. "Peaks at 0.167 vs 0.091" is real but is within the s/B-at-peak spread of Kimura's own dense tests (0.09-0.16) (E2); "too ductile post-peak" should be restated as a late peak (E3).
5. "56 % too soft" is correct as computed; but the model stiffens from the origin (concave-up), an artefact of the p' -> 0 surface with no footing weight (E5, A5, D9). The number is not a clean material verdict.

No error found in: units, s/B, q computation (verified by equilibrium), digitization (all points within 30 kPa), parameter transcription, nu (deliberate nu* device, verified 0.05 in the push).

## Not establishable here; what settles it

- Kimura et al. (1979) (Soils and Foundations / JSCE): roughness of the anisotropy footing, e_max/e_min and G_s of that Toyoura batch, gamma of the poured beds, sand depth, wall lubrication, g reference radius. Zotero MCP failed to connect this session.
- Model side: W_TYR_b16 (running) and a shear:15 leg at B 0.9 for the mesh band; a leg at e_init 0.6586; a leg with 16 kPa footing weight (or the test's actual footing mass); a K0 0.4 leg; an R1-OFF leg that survives past s/B 0.02 (none exists at B 0.9).
