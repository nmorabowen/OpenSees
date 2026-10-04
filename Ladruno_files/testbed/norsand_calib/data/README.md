# WP-144 P3 calibration data pack (Toyoura + Ottawa F65)

Built 2026-10-02 for `Ladruno_implementation/144_ladruno_norsand_plan.md` §3 P3 and §5.2 K3/K5 (re-aimed 2026-10-02).
Everything here is **input data with provenance**: no fitting, no tuning. Derived columns are labelled "derived" in the
CSV headers and in this file. Every number names its source (paper + figure/page, or the file and the build script).

```
data/
  README.md                         this file
  build_data.py                     reproducible builder (reformats + labelled derived columns); run from repo root
  tatsuoka1986/                     K3 calibration data (drained plane strain, Toyoura)
    tats86_fig4a_psl24.csv  tats86_fig6a_d90.csv  tats86_fig16a_iso.csv  tats86_fig8_s3_4p9.csv   one CSV per test
    tests_meta.csv                  per-test metadata + derived peak values
    fig9_phi_peak_vs_e.csv  fig22_eps_peak_vs_e.csv   phi_peak(e), eps_peak(e) series, copied from TIMs
    point_tests.csv                 Fig. 9 + Fig. 22 joined (harness point-test schema)
    source/                         UNEDITED copies of the TIMs files (tests_data.py, tatsuoka_digitized.json, data_f4a.json, cases.json, compare.md)
  kimura1985/                       K5 primary data (centrifuge strip footings, Toyoura)
    source/kimura1985_fig9_digitized.csv   UNEDITED copy of the TIMs digitisation (6 curves, 1404 points)
    source/redteam_report_TIMs_2026-09-29.md   UNEDITED copy of the TIMs red-team of that digitisation (read it: it carries the caveats)
    tests_meta.csv                  per-curve metadata + peaks
  dm04/
    toyoura_table1.csv              DM04 Table 1 (the 15 constants), with equation roles
    toyoura_csl_derived.csv         e_c(p) and psi = e - e_c at the e of the tests (arithmetic, p_at = 100 kPa assumed)
  ottawa_f65/                       K3 methodology check (Ottawa F65, triaxial TC + TE, LEAP-2015)
    ottawa_f65_{TC|TE}_e{724|604|584}_s{kPa}.csv   16 files, columns as in the originals + units header
    tests_meta.csv
    curves/                         (2026-10-03) the same 16 tests as harness lab-curve files (sr from effective stresses, eps_a_sign -1 for TE)
  wang_toyoura/                     (2026-10-03) Wang database, Toyoura HKU (7 tests) + Tokyo (6 tests), harness lab-curve files + tests_meta.csv
  fukushima1984/                    (2026-10-03) Fukushima & Tatsuoka 1984: Table 1 (78 TC tests), Fig. 5 curves, Figs. 17/19/20 summary tables
  lam_tatsuoka1988/                 (2026-10-03) Lam & Tatsuoka 1988: TC / PSC curves (Figs. 5, 9), phi(b, orientation, e) tables (Figs. 4, 8, 21)
  chakraborty2009/                  (2026-10-03) Chakraborty & Salgado 2009: Toyoura TC / PSC Q, R regression tables
  build_lab_curves.py               Wang + Ottawa converter;  build_literature_curves.py  FT/LT converter;  _digitise/  the digitising tool-kit (section 9)
```

## 1. Conventions

- **Units.** Stresses kPa, strains %, area mm^2, load N, time min. Tatsuoka quotes kgf/cm^2; the paper's own conversion is
  1 kgf/cm^2 = 98 kN/m^2 (Fig. 2 caption, p. 67), so 0.05 / 0.1 / 0.5 / 1.0 / 4.0 kgf/cm^2 = 4.9 / 9.8 / 49 / 98 / 392 kPa.
- **Signs, all datasets: compression positive.** Tatsuoka: eps_a > 0 in compression, **eps_v dilation NEGATIVE** (the paper's axis
  convention, kept; CSV key `eps_v_convention: dilation_negative`). Ottawa: `volumetric_strain` and `vertical_strain` compression
  positive, so dilation and axial extension are negative. Kimura: settlement and load positive.
  OpenSees is tension-positive, so a driver must flip the sign.
- **Tatsuoka CSV schema** = `Ladruno_files/testbed/norsand_calib/harness/data.py` lab-curve schema: `eps_a_pct, sr, eps_v_pct`
  with `# key: value` metadata lines (`id, kind, sigma3_kPa, e0, source, eps_v_convention`). The ratio and eps_v curves were
  digitised at different eps_a, so each row fills **one** of `sr` / `eps_v_pct` (the other is empty). Loads with `load_curve`
  (checked 2026-10-02 on Esmeralda: 4 curves + 22 point tests).
- Plane strain: eps_2 = 0, so eps_v = eps_1 + eps_3 (compression positive).
- sr = sigma1'/sigma3' (major over minor effective stress); phi = asin((sr - 1)/(sr + 1)).

## 2. Provenance table

| Dataset | Files | Source (what and where) | How it got here | Modifications | Uncertainty |
|---|---|---|---|---|---|
| Tatsuoka 1986 stress-strain curves, 4 tests | `tatsuoka1986/tats86_*.csv`, `source/tests_data.py`, `source/tatsuoka_digitized.json` | Tatsuoka, Sakamoto, Kawamura, Fukushima (1986) *Soils and Foundations* 26(1):65-84, "Strength and deformation characteristics of sand in plane strain compression at extremely low pressures". Figs. 4(a) p.71, 6(a) p.73, 8 p.73, 16(a) p.78. PDF: TIMs `References/tatsuoka1986_lowpressure_ps.pdf` (journal page = PDF page + 64) | Digitised by the TIMs Workbench (`Tries/2d-model/references/tatsuoka1986_element/`, files dated 2026-09-29). Copied byte-identical (`cmp` checked) into `source/`; reformatted by `build_data.py` | Reformat only (merge ratio and eps_v columns, metadata header). Values unchanged | PSL-24: circle markers auto-detected, 1 px = 0.02 % strain, 0.011 ratio (TIMs docstring). Others: line read by eye on a 0.05-0.1 % grid, resolution not quantified by the digitiser. My visual cross-check (section 5): peak ratio within about 0.05, eps_peak within about 0.1 % |
| Tatsuoka phi_peak(e), eps_peak(e) | `tatsuoka1986/fig9_phi_peak_vs_e.csv`, `fig22_eps_peak_vs_e.csv`, `point_tests.csv` | Figs. 9 (p.75, phi by method C-2-T, delta = 90 deg) and 22 (p.82, eps_a at (sigma1'/sigma3')max, delta = 90 deg), isotropically consolidated, plane strain | `FIG9`, `FIG22` dicts of TIMs `tests_data.py`, copied unedited | Reformat; `point_tests.csv` joins the two by sigma_c' and e within 0.0035 | About 0.3 deg in phi, 0.003 in e, 0.1 % in eps (my visual check). Known gap: sigma_c' = 0.1 series has 2 points in Fig. 9 against more markers on the page (section 6, item 5) |
| Kimura 1985 Fig. 9 | `kimura1985/source/kimura1985_fig9_digitized.csv`, `tests_meta.csv` | Kimura, Kusakabe, Saitoh (1985) *Geotechnique* 35(1):33-45, Fig. 9 (journal p.40 = PDF p.8): strip footing, B = 30 mm, 30 g, Toyoura sand | TIMs `Tries/2d-model/references/kimura1985_fig9/` (digitised 2026-09-29 from a 600-dpi render; the calibration and tracking method are in the CSV header), copied byte-identical | None to the CSV | TIMs: +-0.02 mm settlement (+-0.0007 s/B), +-5 kPa, up to +-0.06 mm where curves overlap (q < 1.5x100 kPa for all, q < 7.5x100 for H86.7/H76.5/V64.6). TIMs red-team: all points within 30 kPa of an independent re-read. Drafting accuracy of the original figure is not included |
| DM04 Toyoura constants | `dm04/toyoura_table1.csv` | Dafalias & Manzari (2004) *J. Eng. Mech.* 130(6):622-634, **Table 1 on journal p. 626 (PDF p.5)**, equations in Table 2 p. 628; PDF TIMs `References/dafalias2004.pdf` | Read visually from a 250-dpi render of the table, then compared with the cluster oracle (section 4) | None | Transcription: none found, all 15 values agree with two independent copies |
| Ottawa F65 monotonic drained triaxial, 16 tests | `ottawa_f65/*.csv` | George Washington University LEAP-2015 database (Vasko, El Ghoraiby, Manzari, Dec 2014): `download.zip` -> `Monotonic Triaxial Experiments/{Loose - eo=0.724, Dense - eo=0.584, Test Density - eo=0.604}/...`. Zip: `C:/Users/nmora/Dropbox/SOILS_rev/download.zip`, 4,665,420 bytes, sha256 108517a3c71cbd57d052f42bb4efd964b3efbdf3d030c7b5c4709f19ac200cb3, never modified | Extracted from the zip by `build_data.py`; the cyclic tests and the PDFs were not copied | The 3-row text header replaced by a `#` header with units; numbers copied as written (8 significant digits), the index column kept | Lab data, no uncertainty given. Characterisation: Gs 2.648 (6 trials), e_max 0.7389 +- 0.0247, e_min 0.4915 +- 0.0183 (9 trials, `Ottawa-F65 Sand Characterization Tests.pdf` slide 7) |
| TIMs red-team of the Kimura digitisation | `kimura1985/source/redteam_report_TIMs_2026-09-29.md` | TIMs Workbench, 2026-09-29 | copied unedited | none | Not re-verified here beyond the checks of section 5 |
| Wang database, Toyoura HKU + Tokyo (13 drained TC) | `wang_toyoura/` | Wang, Sand Triaxial Test Database rev2 (4TU, doi 10.4121/086847a6-ba39-4d66-973b-6b93028c7ad8, CC-BY-4.0), `integrated_dataset/33_Toyoura_sand_HKU.csv` (Chen & Yang 2025, Eng. Geol. 345:107863) and `34_Toyoura_sand_Tokyo.csv` (Verdugo & Ishihara 1996, S&F 36(2):81-91). Zip `C:/Users/nmora/Dropbox/SOILS_rev/WP144_calibration/Wang_sand_triaxial_database_rev2.zip`, 37,381,268 bytes, sha256 4c351c9f5f03a0cf5feff868473c191830e6f42931cbbe489dce8f7c1c134b8b, read only | `build_lab_curves.py` (this round's check: one HKU row compared with the zip by hand, exact) | Reformat only (sr = sigma_v/sigma_r) | The database itself compiled these from published figures (digitised; strains to 4 decimals); tolerance not stated |
| Ottawa F65 as lab curves | `ottawa_f65/curves/` | the 16 files above | `build_lab_curves.py`; sr = max/min of (total - pore) stresses; hand check of TE e604 s100 row 2: 2.357813 | `eps_a_sign` -1 for the 4 TE tests (the file's vertical strain is negative in extension) | none beyond the lab's |
| Fukushima & Tatsuoka 1984 Table 1 | `fukushima1984/table1_tests.csv`, `point_tests_tc.csv` | S&F 24(4):30-48, Table 1 journal p. 32 (PDF page 3), PDF `Fukushima_Tatsuoka_1984_SF24-4_30.pdf` (journal page = PDF page + 29) | transcribed by eye from native-resolution crops of the NII scan, two passes | kgf/cm2 x 98 -> kPa; sr derived | printed numbers: none expected |
| Fukushima & Tatsuoka 1984 Fig. 5 curves | `fukushima1984/curves/` (12 files) | Fig. 5(a), (b), (d), journal p. 36 (PDF page 7) | digitised (section 9) | none | sr +-0.05, eps_a +-0.1 %, eps_v +-0.15 %; check against Table 1 phi: max difference 0.054 in sr |
| Fukushima & Tatsuoka 1984 Figs. 17, 19, 20 | `fukushima1984/fig17_*.csv`, `fig19_*.csv`, `fig20_*.csv` | journal pp. 45-46 (PDF pages 16-17) | Fig. 17 traced, 19/20 markers read | none | phi +-0.1 deg (Fig. 17), v +-0.1 %, eps_a +-0.05 %, sigma3' +-4 % |
| Lam & Tatsuoka 1988 curves | `lam_tatsuoka1988/curves/` (10 files) | S&F 28(1):89-106, Fig. 5 (p. 94, PDF page 6) TC; Fig. 9 (p. 97, PDF page 8) PSC; PDF `Lam_Tatsuoka_1988_SF28-1_89.pdf` (journal page = PDF page + 88) | digitised (section 9) | strain axis converted from eps1-eps3 to eps_a | sr +-0.04-0.05, eps_a +-0.15-0.2 % |
| Lam & Tatsuoka 1988 phi tables | `lam_tatsuoka1988/fig21a_*.csv`, `fig21b_*.csv`, `fig4a_*.csv`, `fig8a_*.csv` | Figs. 21 (p. 101), 4(a) (p. 94), 8(a) (p. 96) | markers read on a gridded overlay | none | phi +-0.25-0.3 deg, e +-0.005, b +-0.02 |
| Chakraborty & Salgado 2009 | `chakraborty2009/` | ICSMGE 17, pp. 652-655, Tables 1-2 (p. 654), Fig. 3 (p. 655); PDF `Chakraborty_Salgado_2009_ISSMGE_STAL0652.pdf` | transcribed from a 250-dpi render (text layer scrambles Table 2) | none | none expected |

## 3. Tatsuoka et al. 1986: the four tests

All: saturated, air-pluviated, fresh Toyoura sand (mean grain size 0.16 mm, U_c 1.46, Gs 2.64, angular to sub-angular; p. 70),
**delta = 90 deg** (bedding plane normal to sigma1, i.e. the air-pluviation direction; the V case of Kimura 1985, p. 39),
isotropic consolidation, drained plane strain compression at 0.25 %/min axial strain rate (p. 70), stresses by method **C-2-T**
(membrane forces and confining-plate friction corrected; the paper adopts this method for all results, pp. 71-72).

| test_id | Fig. / page | sigma3' kPa | e_0.05 | D_r (assumed, 3.2) | R_peak | eps_a at peak % | phi_peak deg (derived) | max dilatancy -d eps_v / d eps_a (derived) | eps_v curve |
|---|---|---|---|---|---|---|---|---|---|
| `tats86_fig4a_psl24` (PSL-24) | 4(a) / 71 | 4.9 | 0.714 | 0.692 | 6.32 | 1.44 | 46.6 | 0.93 | yes (23 pts) |
| `tats86_fig6a_d90` | 6(a) / 73 | 4.9 | 0.700 | 0.729 | 6.52 | 1.74 | 47.2 | 0.94 | yes (11 pts) |
| `tats86_fig8_s3_4p9` | 8 / 73 | 4.9 | 0.755 | 0.584 | 5.76 | 1.99 | 44.8 | n/a | **not digitised** |
| `tats86_fig16a_iso` | 16(a) / 78 | 49.0 | 0.716 | 0.687 | 6.32 | 2.05 | 46.6 | 0.68 | yes (8 pts) |

(R_peak, eps_a at peak, phi and dilatancy are derived from the digitised curves by `build_data.py`; the plan's "sigma3' 4.9 / 49 kPa,
e 0.700-0.755" is exactly this set.)

### 3.1 What e means here (read before using e0)

`e_0.05` is the void ratio at sigma_c' = 0.05 kgf/cm^2 = 4.9 kPa, measured on the thawed sample under vacuum (p. 70; Fig. 9 axis).
For the 49 kPa test (Fig. 16(a)) and for all points at sigma_c' >= 98 kPa the void ratio **at the test confining stress is lower and is
not reported** in the paper. In the pack `e0` always means e_0.05. Paper figure legends: Fig. 4(a) 0.714; Fig. 6(a) delta 90 = 0.700;
Fig. 8 sigma_c' 0.05 = 0.755 (the other four curves of Fig. 8 are 0.746, 0.741, 0.748, 0.752 for 0.1, 0.5, 1.0, 4.0 kgf/cm^2); Fig. 16(a) iso = 0.716.

### 3.2 The D_r assumption (e_max/e_min = 0.977 / 0.597)

**Tatsuoka 1986 does not give e_max or e_min** (I read all 20 pages, none appears). The value 0.977 / 0.597 used by TIMs
(`analyse.py`: `(0.977 - e)/0.38`) is the Verdugo & Ishihara batch. It is recoverable from DM04 itself: DM04 quotes (e, D_r) pairs
0.735 / 63.7 % (Figs. 5 p. 629), 0.833 / 37.9 % (Fig. 6), 0.907 / 18.5 % (Fig. 7; text p. 632) which give
e_max 0.977 and e_min 0.597 (D_r = (0.977 - e)/0.380 reproduces all three). Tatsuoka's Toyoura batch may differ (the TIMs red-team,
D6, cites other batches with e_min 0.597-0.616; **not verified by me**). So `D_r_assumed` is a convenience for comparing with DM04,
not a measured property of Tatsuoka's samples. The state variable that matters for NorSand is psi = e - e_c, not D_r.

### 3.3 "External" strains and digitisation provenance

- eps_a is the **average axial strain from boundary (external) displacements** (Fig. 22 and text p. 82; the stress is from loads at the
  boundaries). The paper says explicitly (p. 74) that the strain field in the sample is not uniform, that shear bands start to form
  before the peak, and that the non-uniformity is larger at low sigma_c'; so the post-peak and the strain at peak are not pure
  material properties of an element. Use the pre-peak branch and the peak value to calibrate; treat the post-peak softening as
  a specimen-with-shear-band response (this is the K3 sanity gate, not a fit target).
- The first point (0, 1) of each CSV is a TIMs prepend. In Fig. 6(a) the sigma1'/sigma3' curves leave the axis well above 1 (about 2 for
  the dense delta = 90 curve, my reading) and rise steeply within about 0.1 % strain; p. 79 attributes the sigma2'/sigma3' > 1 at the start of
  shearing to seating stress introduced when setting the confining plates. Do not fit E_sec below 0.1 % without reading the figure.
- PSL-24 has several rows with the same eps_a and different R in the pre-peak part (circle markers overlapping, 1 px = 0.02 %). The
  TIMs code made the abscissa non-decreasing before using it. The file keeps the raw points; `load_curve` sorts stably by eps_a.
- Digitised by TIMs 2026-09-29 (files in `source/`; `data_f4a.json` is the raw auto-detection of Fig. 4(a), `compare.md` is TIMs'
  DM04-versus-Tatsuoka table, kept because it states the DM04 plane-strain misfit that P3 starts from: DM04 gives eps_peak 0.46-0.58x of
  the test at sigma3' 4.9 kPa and about 1.0x at 49 kPa).

### 3.4 Series for phi_peak(e), eps_peak(e)

`fig9_phi_peak_vs_e.csv` (22 points, delta = 90 deg, 5 stress levels) and `fig22_eps_peak_vs_e.csv` (22 points) are the TIMs
`FIG9` / `FIG22` dictionaries, converted to kPa. Ranges: phi_peak 37.5-49.9 deg at e_0.05 0.862-0.653; eps_peak 1.14-5.43 %.
`point_tests.csv` carries the pairs (phi, eps_peak) for the same sample where the two figures agree in e within 0.0035; the other points
keep an empty eps_peak.

## 4. DM04 Toyoura Table 1 and comparison with the cluster oracle

| Group | Constant | Value | Role (DM04 Table 2, p. 628) |
|---|---|---|---|
| Elasticity | G0 | 125 | G = G0 p_at (2.97 - e)^2/(1 + e) (p/p_at)^0.5 |
| Elasticity | nu | 0.05 | K = 2(1+nu)G/(3(1-2nu)) |
| Critical state | M | 1.25 | CSL slope in triaxial compression |
| Critical state | c | 0.712 | M_e/M_c |
| Critical state | lambda_c | 0.019 | e_c = e0 - lambda_c (p_c/p_at)^xi |
| Critical state | e0 | 0.934 | CSL intercept (not the initial e) |
| Critical state | xi | 0.7 | CSL exponent |
| Yield surface | m | 0.01 | yield-cone opening |
| Plastic modulus | h0 | 7.05 | b0 = G0 h0 (1 - c_h e)(p/p_at)^-0.5 |
| Plastic modulus | c_h | 0.968 | idem |
| Plastic modulus | n^b | 1.1 | M^b = M exp(-n^b psi) |
| Dilatancy | A0 | 0.704 | A_d = A0 (1 + <z:n>) |
| Dilatancy | n^d | 3.5 | M^d = M exp(n^d psi) |
| Fabric-dilatancy | z_max | 4 | cyclic only |
| Fabric-dilatancy | c_z | 600 | cyclic only |

Source: DM04 Table 1, journal p. 626 (PDF page 5). Calibrated by DM04 to Verdugo & Ishihara (1996) triaxial tests, **p' 100-3000 kPa,
D_r 18.5-63.7 %, e 0.735-0.907** (text p. 632). **p_at is not tabulated.**

**Comparison with `gate0_toyoura_oracle.py`** (TIMs `Tries/2d-model/model/C4-sanisand-fork/regularization-wp150/`, read only):

- Its docstring (lines 3-5) lists G0 125, nu 0.05, M 1.25, c 0.712, lambda_c 0.019, e0 0.934, xi 0.7, m 0.01, h0 7.05, ch 0.968, nb 1.1, A0 0.704,
  nd 3.5, zmax 4, cz 600, claiming to be "verified against the PDF, p. 5 = journal p. 626". **All 15 agree with Table 1.**
- The code constants it imports (`TOYOURA`, `sanisand_reference/model.py` lines 147-149, identical in this worktree's
  `Ladruno_scripts/sanisand_reference/`, in TIMs `oracle-wp134/` and in `reseat-r1-wp151/`) are the same vector
  `[125, 0.05, 0.8, 1.25, 0.712, 0.019, 0.934, 0.7, 100, 0.01, 7.05, 0.968, 1.1, 0.704, 3.5, 4.0, 600]` (argument order of
  `LadrunoSANISAND`: G0 nu e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz).
- **No difference in the DM04 constants.** The only items that are not in Table 1: **P_atm = 100** (an explicit choice in the oracle;
  the TIMs campaign set uses 101), **e_init** (0.8 is a placeholder in `TOYOURA`; gate0 overrides it per test with the test's e0), and
  `Den` 2.0 (mass density, not constitutive). The paper's formulation (`paper` option) and the UW options (`uw_model`: G from e_init,
  low-p D sigmoid, p_min) differ as documented in the oracle; Table 1 is the same for both.
- The TIMs campaign set (G0 264.32, M 1.3309, lambda_c 0.027, e0 0.83, xi 0.45, nb 3.5, A0 0.05, nd 5.75, ...) is a different sand, not Toyoura
  (plan §3 P3). Not used here.

Derived arithmetic (script, p_at = 100 kPa assumed; `dm04/toyoura_csl_derived.csv`):

- phi_c = asin(3M/(6 + M)) = 31.15 deg from M = 1.25; M_e = c M = 0.890 gives phi_e = asin(3 M_e/(6 - M_e)) = 31.50 deg (equal friction in
  compression and extension would need c = 0.706). Not a defect, only the size of the DM04 rounding.
- e_c(p) = 0.934 - 0.019 (p/100)^0.7: 0.9317 at 4.9 kPa, 0.9225 at 49 kPa, 0.915 at 100 kPa, 0.8754 at 500 kPa.
  psi at Tatsuoka's e_0.05 (0.700-0.755) is -0.23 to -0.18 at 4.9 kPa and -0.22 to -0.17 at 49 kPa (the void ratio at 49 kPa is lower than e_0.05,
  so psi there is somewhat more negative). Kimura V85.6 (e about 0.652 under the 3.2 assumption) is psi = -0.28 at 4.9 kPa.
  **The power-law CSL is an extrapolation by one to two decades of p below the 100-3000 kPa range it was calibrated on**; whether
  e_c(5 kPa) = 0.93 is right is not testable with the Tatsuoka data (loose samples e_0.05 0.80-0.86 still dilate at 4.9 kPa, p. 73 Fig. 6(b), which only says e_c(4.9 kPa) is above 0.80).
- G at e = 0.70, p = p_at: G = 125 x 100 x (2.97 - 0.70)^2/(1 + 0.70) = 37.9 MPa (p_at = 100; 101.3 changes this by 0.65 %).
- The fork's NorSand CSL takes `-p_a`; use 100 kPa for continuity with the oracle.

## 5. Cross-check against the PDFs (done 2026-10-02, visual, 250-dpi renders)

Rendered with PyMuPDF on Esmeralda (`~/wp144/data_build/`, scratch) and read by eye; values are my readings, not digitisations.

| Item | Read from the PDF | Pack value | Verdict |
|---|---|---|---|
| Fig. 4(a), p.71, PSL-24 | peak R about 6.35 at eps_a about 1.5 %; R about 4.33 at 5 %, about 4.5 at 10 %; eps_v about -2.7 % at 10 % (axis -2 / -4 / -6) | 6.32 at 1.44 %; 4.31-4.36 at 5.3 %; 4.47 at 9.5 %; eps_v -2.76 at 9.5 % | agrees (R within 0.05, eps within 0.1 %, eps_v within 0.1) |
| Fig. 6(a), p.73, delta 90 | peak R about 6.54 at about 1.7 %; line ends R about 4.7 at 8.8 %; legend e 0.700 | 6.52 at 1.74 %; 4.68 at 8.83 %; 0.700 | agrees |
| Fig. 8, p.73 | sigma_c' 0.05 curve peak about 5.8 at about 2 %; legend e 0.755 | 5.76 at 1.99 %; 0.755 | agrees |
| Fig. 16(a), p.78 | iso curve peak about 6.3 at about 2 %; legend e 0.716 | 6.32 at 2.05 %; 0.716 | agrees |
| Fig. 9, p.75 | every plotted marker of the 0.05, 0.5, 1.0, 4.0 series | `fig9_phi_peak_vs_e.csv` | agrees to about 0.3 deg / 0.003 in e. The 0.1 series is incomplete (section 6.5) |
| Fig. 22, p.82 | 0.05 series (0.755 -> 2.05, 0.715 -> 1.45, 0.700 -> 2.0), 4.0 series (0.674 -> 3.62, 0.752 -> 4.45, 0.82 -> 5.45) | `fig22_eps_peak_vs_e.csv` | agrees to about 0.1 % |
| Kimura Fig. 9, p.40 | peaks (x100 kN/m^2 at S mm): V85.6 19.5 at 2.7; H86.7 17.3 at 4.9; V75.1 12.4 at about 3.2; H76.5 11.6 at about 4.0; V64.6 8.1 at 2.5; H61.0 7.9 at about 6.1 | 1953 at 2.72; 1734 at 4.89; 1236 at 3.14; 1161 at 3.92; 814 at 2.46; 792 at 6.01 (kPa at mm) | agrees (within 1 %, S within 0.1 mm) |
| Kimura Table 3, p.42 (independent) | observed N_gamma = 174 for the V case; 174 x 0.5 x 15.9 x 0.9 = 1245 kPa | V75.1 peak 1236 kPa (-0.7 %) | agrees (TIMs red-team A1 found the same) |
| Ottawa monotonic PDF, slides 4-11 | q-eps_a and eps_v-eps_a for loose / dense / test density, TC and TE | files | shape, peaks and signs agree (e.g. dense TC 100 kPa q peak about 390 at 5 %: file 390.6 at 4.7 %) |

Not checked (stated, not hidden): the individual points of Fig. 6(a), Fig. 8 and Fig. 16(a) below the peak; I checked only peak, end point and legend.

## 6. Inconsistencies and caveats register

1. **e_max/e_min is not in Tatsuoka 1986** (section 3.2). D_r is assumed; use e and psi.
2. **e is e_0.05, not the void ratio at the test stress** (section 3.1). For the 49 kPa test the model e0 should be e_0.05 minus a compression
   increment that the paper does not give; the plan's K3 range (e 0.700-0.755) is e_0.05.
3. **Fig. 8 eps_v curve, and the other Fig. 8 curves (sigma_c' 0.1-4.0 kgf/cm^2, e about 0.74-0.75), are not digitised.**
   Fig. 7(a) (sigma_c' 4.0, dense, delta 90 e 0.714) and the other delta curves are also not digitised. They would give the h(p) trend
   across 4.9-392 kPa at one density; a candidate for a follow-up digitisation, not done here (out of scope, read-only on TIMs).
4. Digitisation resolution of the three line-read tests is not quantified by TIMs; PSL-24 has 1-px duplicates (section 3.3).
5. **Fig. 9 sigma_c' = 0.1 kgf/cm^2 series has 2 points in `FIG9`**, while Fig. 22 gives 4 for the same stress level and the Fig. 9 page shows further "x"
   markers (near e 0.746 and 0.805, phi about 44 and 41 deg) that were not digitised. Fig. 9 and Fig. 22 e values for the same sample differ by up to
   about 0.005 (e.g. 0.672 vs 0.677 for the densest 4.9 kPa tests), so `point_tests.csv` leaves one 4.9 kPa eps_peak empty rather than guess.
6. **Kimura Fig. 9 (from the TIMs red-team, retained as caveats; I verified items 6a-6c against the paper text):**
   a. The paper does not give e_max, e_min, Gs or the unit weight for the pouring series (anisotropy series, p. 39). e and gamma in
      `kimura1985/tests_meta.csv` are **assumed**: e = 0.977 - D_r x 0.38 (V85.6: 0.6517, the value TIMs used), gamma = 2.64 x 9.81/(1 + e) (V85.6: 15.68 kN/m^3).
      The TIMs deck used gamma = 15.9 kN/m^3, which is the vibro-compacted de Beer series of p. 36, about 1 % too high for e = 0.6517.
   b. The footing roughness of the anisotropy series is not stated (roughness is stated for the roughness and de Beer series, pp. 36-37). On p. 37 the
      smooth footing gives 80 % of the rough capacity at B_N 1.2 m. If the Fig. 9 footing was smooth, a rough-footing model is 25-35 % high.
   c. Which centrifuge produced Fig. 9 is not stated. Mark I container (Table 1, p. 34): 0.5 m long, 0.1 m wide, 0.3 m deep, effective radius 1.18 m;
      Mark II: effective radius 1.25 m (Table 2). A 0.1 m wide box (L/B = 3.3) is not perfectly plane strain; side-wall friction typically raises the
      measured capacity (the red-team estimates 5-15 %, unverified).
   d. The 30 g field is radial (about +-13 % over the 0.3 m depth); embedment is not stated for Fig. 9 (assumed surface).
   e. D_r of the V and H curves are the paper's labels (to 0.1 %); the curve-to-label assignment of V75.1, H76.5, V64.6, H61.0 rests on the leader lines of the figure
      and on my reading of the page (section 5); TIMs verified the V85.6 and H86.7 assignments independently, and V75.1 through the Table 3 N_gamma.
   f. Kimura's Toyoura (plane strain phi about 49 deg, p. 41) and Tatsuoka's delta = 90 deg samples (phi_peak 44.7-49.9 deg, Fig. 9) are consistent
      with each other; the V case is the delta = 90 deg case.
7. **DM04 versus the cluster oracle: no difference in the 15 constants** (section 4). Items outside Table 1 (P_atm, e_init, Den) are choices, listed there.
8. **Ottawa F65 (this pack):**
   a. The loose tests (e0 0.724) have D_r 0.06 with the averaged e_max/e_min (0.7389 +- 0.0247): at the level of the e_max scatter, i.e. the loosest reproducible state.
   b. Dry-density labels vs e: rho_d 1537 / 1673 / 1652 kg/m^3 reproduce e 0.724 / 0.584 / 0.604 with Gs = 2.65 (e = 2.65/rho_d - 1: 0.7242, 0.5840, 0.6041), not with the
      measured Gs = 2.648 (0.7229, 0.5829, 0.6029). The e labels are used as given (the files and folder names use them).
   c. **TE e0 0.604, sigma3' 100 kPa looks anomalous**: derived phi_peak 44.1 deg against 35.0 (200 kPa) and 36.4 (300 kPa) for the same density and TC 37.4-38.7;
      q is only about -80 kPa (the measurement is at the low end of the load cell), and the file stops at eps_a = -8.1 % (82 rows). Treat it as low-weight or exclude it from a rho fit.
   d. The loose TC tests at 500-700 kPa have not reached a plateau (peak at or near the end of the test, 700 kPa stops at 13.1 %); phi_peak 35.1-35.7 deg at those stresses is a
      lower bound of the sample capacity.
   e. Only one TE test at e0 0.724 (200 kPa); TE tests exist at three stresses only for e0 0.604, which is the pair for the rho-from-data check (TC 100/200/300 and TE 100/200/300).
   f. Triaxial TC/TE, not plane strain: used as a methodology check, not as calibration data for the Toyoura-based strip (plan K3).
   g. The slide "Test Density" extension volumetric plot writes rho_d in g/cm^3 (a typo of the source; the file name and the other slides use kg/m^3).
   h. Total stresses in the files: effective = total - pore pressure (back pressure about 495-555 kPa). Row-1 sigma3' = horizontal - pore reproduces the nominal 100-700 kPa within 0.15 kPa (meta column `sigma3_eff_start_kPa`).
9. **No undrained or extension lab data for Toyoura** exist in this pack (plan K3 states the limit). Ottawa is the only extension data and is a different sand.

Derived Ottawa values (from `ottawa_f65/tests_meta.csv`; phi_peak from effective stresses, s1' the larger of vertical/horizontal):

| e0 | mode | sigma3' kPa : phi_peak deg (eps_a at peak %) |
|---|---|---|
| 0.724 (loose) | TC | 100: 37.3 (16.3); 200: 35.9 (13.5); 300: 35.0 (13.2); 500: 35.1 (19.8); 600: 35.3 (16.0); 700: 35.7 (12.8, test ends 13.1) |
| 0.724 | TE | 200: 32.4 (-6.5) |
| 0.604 (test density) | TC | 100: 38.7 (4.8); 200: 37.5 (5.0); 300: 37.4 (8.6) |
| 0.604 | TE | 100: 44.1 (-3.6, anomalous); 200: 35.0 (-3.3); 300: 36.4 (-3.7) |
| 0.584 (dense) | TC | 100: 40.8 (5.3); 200: 40.0 (5.6); 300: 39.9 (5.6) |

## 7. Upstream notes (for TIMs; not edited here)

**TIMs e_init inconsistency.** The initial void ratio used for the SANISAND campaign sand is **0.6944** in `Tries/2d-model/model/strip.py` line 154
(`SANISAND_COMMON["e_init"]`; the same value is in the cluster's footing deck `footing_ab.py:90`, `sanisand_reference` `CAMPAIGN`, and `regularization-wp150/mono_fit.py`)
but **0.704** in `labs/sfim/SOIL.py` line 36 (`MD04['e0'] = 0.704`, comment "Initial void ratio"). Every other MD04 constant listed in `SOIL.py` (G0 264.32,
c 0.71, lambda_c 0.027, ec0 0.83, ksi 0.45, P_atm 101, m 0.005, h0 1.3, ch 0.968, nb 3.5, A0 0.05, nd 5.75, z_max 12.5, cz 1100; nu = K0/(1+K0) from phi = 33 deg = 0.3129;
Mc from phi = 33 deg = 1.3309) equals the `strip.py` value, so the 0.010 difference is isolated. Effect: e_init enters only through G (UW option: G with e_init;
(2.97 - e)^2/(1 + e) is 3.056 at 0.6944 and 3.013 at 0.704, 1.4 % lower). It is an upstream inconsistency to settle on the TIMs side; this pack does not use either value.

## 8. Rebuild

```
python Ladruno_files/testbed/norsand_calib/data/build_data.py <path to download.zip>
```
(run on Esmeralda per the WP-144 rule; needs numpy). It reads `source/` and the zip, writes every CSV above except the `source/` copies, and prints the DM04 phi_c / phi_e.
Cross-checks of section 5 were done by eye and are not scripted.

## 9. P3 data additions, 2026-10-03 (WP-144 literature round)

Everything below follows the conventions of section 1 (stresses kPa, strains %, compression positive, `eps_v` dilation negative, sr = sigma1'/sigma3'). No fitting, no tuning;
"derived" columns are labelled. Curve files load with `harness/data.py:load_curve` (checked 2026-10-03: 22 new literature curves, 13 Wang and 16 Ottawa curves, 78 point tests).

### 9.1 What was added

| Dataset | Files | Content |
|---|---|---|
| Wang Toyoura (HKU, Tokyo) | `wang_toyoura/wang{33_HKU,34_Tokyo}_TMD<n>.csv`, `tests_meta.csv` | 13 drained isotropically consolidated TC tests, written by `build_lab_curves.py` |
| Ottawa F65 lab curves | `ottawa_f65/curves/*.csv`, `curves_meta.csv` | 16 tests (12 TC, 4 TE), sr from effective stresses, TE `eps_a_sign = -1` |
| Fukushima & Tatsuoka 1984 | `fukushima1984/table1_tests.csv`, `point_tests_tc.csv` | all 78 TC tests of Table 1: sigma_c', (sigma3')_f, e_0.3, t0, phi (uncorrected); `point_tests_tc.csv` = harness point-test schema |
| | `fukushima1984/curves/ft84_fig5{a,b,d}_sc*.csv` (12) | stress-ratio and eps_v curves of Fig. 5(a) (dense, sigma_c' 0.1-4.0 kgf/cm2: 6), (b) (dense, 0.02-0.1: 3), (d) (loose, 0.02-0.1: 3) |
| | `fig17_phi_corrected_vs_sigma3.csv`, `fig19_*`, `fig20_*` | corrected phi(sigma3') at e_0.3 0.70 / 0.85 (2 curves, 44 + 44 points); eps_v at eps_a 5 and 10 % vs sigma3' (31 markers); eps_a at sigma1'/sigma3' = 3, 4 vs sigma3' (24 markers) |
| Lam & Tatsuoka 1988 | `lam_tatsuoka1988/curves/lt88_fig5{a,b}_*` (8), `lt88_fig9_psc_*` (2) | TC (H/W 1.0 and 0.25, omega 0/30/60/90) and PSC (omega 0, H/W 1.9 and 1.0) stress-ratio / eps_v curves |
| | `fig21a_phi_vs_orientation.csv`, `fig21b_phi_vs_b.csv`, `fig4a_tc_phi_vs_e.csv`, `fig8a_psc_phi_vs_e.csv` | phi vs bedding orientation (TC, PSC, TE), phi vs b = 0 / 0.25-0.36 / 1 (TC / PSC / TE), phi(e) for TC and PSC |
| Chakraborty & Salgado 2009 | `chakraborty2009/table1_tx_Q_R.csv` (7 rows), `table2_psc_Q_R.csv` (5), `fit_parameters.csv` | Q, R regression tables and the phi_c / A_psi / Q(sigma_c') constants |

### 9.2 Ranges

- **Wang HKU (7 tests, Chen & Yang 2025):** sigma3 = 100, 200 (x4), 400 (x2) kPa; e0 0.640-0.834 (D_r 38.4-90.6 %); eps_a 0-12 % (121 points, 0.1 % step); phi_peak 31.9-41.4 deg (derived); d50 0.22 mm, Gs 2.64, e_max 0.977, e_min 0.605, phi_c 31.4 (database values).
- **Wang Tokyo (6 tests, Verdugo & Ishihara 1996):** sigma3 = 100 (x3), 500 (x3) kPa; e0 0.810-0.996 (D_r 44 to -5 %, TMD6 has e0 above the reported e_max); eps_a to 26-32 % (262-317 points); phi_peak 30.6-33.1 deg; two loose tests (TMD3, TMD6) have no interior peak. d50 0.17 mm, Gs 2.65, e_max 0.977, e_min 0.597, phi_c 31.
  **No Wang test has sigma3 below 100 kPa**; the low-stress end of the pack is Fukushima & Tatsuoka (sigma_c' 1.96-392 kPa).
  Both files say `drained`, `triaxial compression`, `isotropic` in their own TEST_PROGRAMME; stresses are used as given (sr = sigma_v/sigma_r). The HKU curves come from a small-strain stiffness study (sr is already 1.5 at 0.1 % strain).
- **Fukushima & Tatsuoka Table 1:** (sigma3')_f 2.55-197 kPa measured (sigma_c' 1.96-392 kPa; the 392 kPa tests have no measured value), e_0.3 0.650-0.908, phi (uncorrected) 32.9-47.2 deg. **e_0.3 is the void ratio at 29.4 kPa, not at the shear stress** (paper p. 35, Fig. 4): for the other sigma_c' the shear void ratio differs and is not reported.
- **Lam & Tatsuoka:** all at sigma3' = 98 kPa, e_0.3 0.637-0.675 for the curves (0.70 nominal for the phi tables); H/W 1.9 (PSC) and 1.0 (TC) are the nominal specimens.

### 9.3 How the figures were digitised (tool-kit in `_digitise/`)

1. The two Soils and Foundations PDFs are NII scans: one embedded 1-bit bitmap per page, 2848 px wide (about 345 dpi). `extract_native_images.py` pulls these bitmaps out unchanged (no resampling).
2. **Axis calibration per panel:** the frame lines are fitted (robust line fit through ink runs; the scans are skewed by 0.3-0.9 deg and the skew differs between the left and right frame lines, so a pure rotation is not enough), the tick marks along the
   frame are detected, and data coordinates are mapped by bilinear interpolation between the left and right frame lines (x: values at the two frame lines, y: linear fit through the tick rows of each line; residuals 0.2-3 px, `Panel.resL/resR`). Log axes use log10 as the abscissa.
3. **Curve extraction:** marker blobs are found after filling small holes (open markers) and a morphological opening; curves are followed (a) by a marker chain in pixel space, (b) by a column-cluster family tracker in which every ink cluster can be claimed by one curve only, (c) by a row-wise tracker for steep starts (Fig. 5(d) lower curve), (d) by a single-curve follower for isolated lines (Fig. 17). Seeds are marker positions read on gridded zooms. Every result was inspected on an overlay of the extracted points on the scan, outliers (legend text, neighbouring curves) were removed (`curate.py` despike, x-ranges) and the surviving points are the files. Where markers of several curves overlap into one black band the curve is **not** digitised there (section 9.5); nothing is interpolated into the files except the strain-axis conversion of Lam & Tatsuoka (below).
4. Fig. 9 (PSC): the six series with many similar markers were assigned marker by marker (blob IDs in `_digitise/run_lt9b.py`) after viewing numbered overlays.
5. Table 1, the Chakraborty tables and the b / orientation / e tables are read directly (Fig. 21 and 4/8 markers on a gridded zoom).

### 9.4 Checks and uncertainty

- **Fukushima & Tatsuoka Fig. 5 against Table 1 (independent):** the digitised maximum sr of each curve was compared with (1 + sin phi)/(1 - sin phi) of the Table 1 test with the same sigma_c' and e_0.3 (legend). Differences (12 curves): +0.036, -0.015, +0.011, -0.004, +0.005, -0.018 (5a); -0.031, -0.026, -0.042 (5b); +0.031, -0.053, +0.002 (5d); all within 1.3 % of sr. This validates the sr axis calibration to about +-0.05. The same check is stored per curve in `curves_meta.csv` (`dsr`).
- **Lam & Tatsuoka internal checks:** PSC H/W 1.9 peak sr 6.494 (phi 47.15 deg) against Fig. 8(a) phi at e 0.657: 47.16 deg; H/W 1.0: 47.68 against 47.43 deg (0.25 deg); TC omega = 0, e 0.671: peak sr 5.494 (43.79 deg) against Fig. 4(a) about 43.6 deg at e 0.669; TC / PSC phi at omega = 0 from Fig. 21(b) markers (41.9 / 45.4 deg) against Fig. 21(a) (41.7 / 45.3 deg).
- **Fig. 17 against the paper's text:** the extrapolations to sigma3' = 0 are quoted as 42.4 and 36.1 deg; the traced curves start at 42.2 and 36.1 deg.
- **Per-point uncertainty (1 sigma, my estimate from pixel size, marker size and the checks above):** FT Fig. 5: sr +-0.05 (+-0.1 where neighbouring markers merge), eps_a +-0.1 %, eps_v +-0.15 %; LT Fig. 5: sr +-0.04, D = eps1-eps3 +-0.1 % (eps_a +-0.15 % after conversion), eps_v +-0.1 %; LT Fig. 9: sr +-0.05, D +-0.15 % (eps_a +-0.2 %), eps_v band +-0.15 %; phi tables +-0.25-0.3 deg, b +-0.02 (PSC) / +-0.01, e +-0.005; Fig. 19/20: v +-0.1 %, eps_a +-0.05 %, sigma3' +-4 %. They are in the `uncertainty` header line of every curve file and do not include the experimenters' own error (membrane forces, bedding error, end restraint).
- **Lam & Tatsuoka strain axis:** the paper plots against eps1 - eps3, not eps_a. `build_literature_curves.py` converts with eps_a = (2D + eps_v)/3 (TC) and (D + eps_v)/2 (PSC) using the same test's eps_v trace (documented in each file header). An error of 0.1 % in eps_v moves eps_a by 0.03-0.05 %.
- **Fukushima & Tatsuoka membrane forces:** Fig. 5 stress ratios are the mid-height values UNCORRECTED for membrane forces (t0 = 0.3 mm in all of 5(a), (b), (d)); at sigma_c' <= 0.1 kgf/cm2 the paper shows this biases phi (Figs. 7, 12-17), and at sigma_c' 0.02-0.05 the actual sigma3' at failure was up to twice sigma_c' (`sigma3_f_kPa` in each header). Use `sigma3_f_kPa` for any p-dependent comparison of those tests.

### 9.5 Not digitised / not readable (stated, not hidden)

1. **FT Fig. 5(c) (loose, sigma_c' 0.1-4.0 kgf/cm2):** the stress-ratio curves of 0.5-4.0 differ by less than 0.1 and their markers merge into one black band; the eps_v curves of the bunch are not separable either. Not digitised. Peak sr and phi of those six tests are in Table 1 (tests 32, 42, 50, 62, 71, 77).
2. **FT Fig. 5(a):** the rise of sr below eps_a about 0.7-1.4 % (the six curves overlap) and the converging tail beyond about 12 % are not given (first / last rows in `curves_meta.csv`); eps_v of sigma_c' 0.1 and 0.2 only up to 7 % (the three upper eps_v curves merge above that; the 0.5 curve continues to 13.5 %).
   **FT Fig. 5(b), (d):** the eps_v traces of sigma_c' 0.02 and 0.05 in 5(b) are one thick pair (same content in both files, flagged in the header); in 5(d) the 0.05 and 0.1 tests coincide below about 9.5 % (the 0.1 file continues alone to 15 %, with real dips at 11.8 and 14.7 % that are in the figure).
3. **LT Fig. 5(a), (b):** the initial rise of sr (D below about 2 %) is not resolved (curves overlap); in 5(b) the omega = 60 and 90 stress-ratio curves coincide above D about 5 % and one merged trace is written to both files.
4. **LT Fig. 9 (PSC):** only the triangle (H/W 1.9) and circle (H/W 1.0) series; sigma1'/sigma3' from D about 1.4 % (15 markers each); sigma2'/sigma3' only 4 markers each (thick overlapped band for D = 6-11 %); eps_v is the centre line of the overlapped band for D below 8 %. The square (0.5) and cross (0.25) series are not digitised (end restraint).
5. **LT triaxial extension (Figs. 12-20):** no TE stress-strain curve was digitised (panels share frames and legends, failure modes make them specimen-dependent); the TE phi at b = 1 is in `fig21*` (single-intersection mode and the Fig. 16 values). Figs. 13, 14, 18 (omega / xi variants of PSC, strains at failure) are not digitised.
6. **FT:** Figs. 6-9 and 12-16 (Table 1 phi plotted against e, membrane-correction comparisons) are not repeated; the membrane-corrected phi per test is NOT tabulated in the paper (only the averaged curves of Fig. 17 and the plotted points of Figs. 14-16, which were not digitised). Fig. 18 (eps_v vs e at eps_a 5 / 10 %) not digitised.
7. **Void ratio at the test stress** is not reported by either paper (e_0.3 only); membrane thickness (FT) and bedding error (LT) are not corrected in any file.
8. The other PDFs in `WP144_calibration/` (Yoshimine & Ishihara 1998, Verdugo & Ishihara 1996 beyond the Wang copy, Siddiquee, Andrade & Ellison, Castonguay & Konrad) were not part of this round.

### 9.6 Rebuild

```
python Ladruno_files/testbed/norsand_calib/data/build_lab_curves.py <Wang_sand_triaxial_database_rev2.zip>   # Wang + Ottawa lab curves
python Ladruno_files/testbed/norsand_calib/data/fukushima1984/build_table1.py
python Ladruno_files/testbed/norsand_calib/data/fukushima1984/build_summary_figs.py
python Ladruno_files/testbed/norsand_calib/data/lam_tatsuoka1988/fig21_phi_tables.py
python Ladruno_files/testbed/norsand_calib/data/lam_tatsuoka1988/fig4_fig8_phi_vs_e.py
python Ladruno_files/testbed/norsand_calib/data/chakraborty2009/build_tables.py
python Ladruno_files/testbed/norsand_calib/data/build_literature_curves.py     # FT/LT curves from the digitised-point JSON files
```
The digitisation itself (`_digitise/`) needs PyMuPDF, scipy and Pillow and the two PDFs (`extract_native_images.py` first); its raw outputs are in `_digitise/outputs/` and, after curation, in the two `digitised_points_*.json` files that the builder reads.
Page images are not committed (copyright of the Japanese Geotechnical Society); only numbers.
