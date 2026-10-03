"""WP-144 P3 data pack builder (reproducible; run from the repo root).

    python Ladruno_files/testbed/norsand_calib/data/build_data.py <path to download.zip>

Reads   data/tatsuoka1986/source/tatsuoka_digitized.json   (copy of TIMs' digitisation, unedited)
        data/kimura1985/source/kimura1985_fig9_digitized.csv (copy, unedited)
        the Ottawa-F65 LEAP-2015 zip (read only; never modified)
Writes  the normalised CSVs next to them (see README.md for the provenance table).

Nothing here fits or tunes anything: it only reformats, and computes the few derived columns that the README
labels as derived (peak values, phi_peak from the ratio, D_r from stated e_max/e_min, the DM04 CSL arithmetic).
"""
import csv
import io
import json
import math
import os
import re
import sys
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ZIP = sys.argv[1]
E_MAX, E_MIN = 0.977, 0.597  # Toyoura, Verdugo & Ishihara batch; implied by the (e, D_r) pairs of DM04 Figs. 5-7
KGF = 98.0  # Tatsuoka et al. 1986 Fig. 2 caption: 1 kgf/cm^2 = 98 kN/m^2


def w(path, header_lines, cols, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        for h in header_lines:
            f.write("# " + h + "\n")
        cw = csv.writer(f, lineterminator="\n")
        cw.writerow(cols)
        for r in rows:
            cw.writerow(r)


def fmt(x, n=6):
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{n}g}"


# ---------------------------------------------------------------- Tatsuoka 1986
D = json.load(open(os.path.join(HERE, "tatsuoka1986", "source", "tatsuoka_digitized.json")))
IDS = {"PSL-24 (Fig. 4a)": "tats86_fig4a_psl24", "Fig. 6(a) delta90": "tats86_fig6a_d90",
       "Fig. 16(a) iso": "tats86_fig16a_iso", "Fig. 8 s3=4.9": "tats86_fig8_s3_4p9"}
NOTE = {"PSL-24 (Fig. 4a)": "circle markers of the original figure, auto-detected (TIMs); +-1 px = 0.02 % strain / 0.011 ratio; "
        "ratio abscissae are not strictly increasing in the pre-peak part because of overlapping circles",
        "Fig. 6(a) delta90": "solid line read by eye on a 0.05-0.1 % grid (TIMs); resolution not quantified by the digitiser",
        "Fig. 16(a) iso": "solid line read by eye on a 0.05-0.1 % grid (TIMs); resolution not quantified by the digitiser",
        "Fig. 8 s3=4.9": "solid line read by eye on a 0.05-0.1 % grid (TIMs); eps_v curve NOT digitised"}
meta = []
for key, tid in IDS.items():
    t = D["TESTS"][key]
    ra = np.array(t["ratio"], float)
    ev = np.array(t["ev"], float) if t["ev"] else np.zeros((0, 2))
    # one CSV per test in the harness lab-curve schema (harness/data.py): eps_a_pct, sr, eps_v_pct; the two curves were
    # digitised at different eps_a, so each row fills one of sr / eps_v_pct and leaves the other empty.
    hdr = [("id", tid), ("kind", "PS"), ("sigma3_kPa", t["s3"]), ("e0", t["e"]),
           ("e0_definition", "e_0.05 = void ratio at sigma_c' = 0.05 kgf/cm2 = 4.9 kPa; the void ratio at the test sigma_3' is NOT reported"),
           ("source", f"Tatsuoka, Sakamoto, Kawamura & Fukushima (1986) Soils and Foundations 26(1):65-84, Fig. {t['fig']}, journal p. {t['page']}"),
           ("eps_v_convention", "dilation_negative"),
           ("material", "saturated air-pluviated Toyoura sand, isotropic consolidation, delta = 90 deg, drained plane strain compression (eps_2 = 0)"),
           ("D_r_assumed", round((E_MAX - t["e"]) / (E_MAX - E_MIN), 3)),
           ("D_r_basis", "(0.977 - e_0.05)/(0.977 - 0.597); e_max/e_min are not in the paper (README)"),
           ("stress_method", "C-2-T (membrane and plate-friction corrected)"),
           ("strain_type", "eps_a = average axial strain from EXTERNAL boundary displacements; eps_v from specimen volume change"),
           ("digitisation", "TIMs Workbench tests_data.py; " + NOTE[key]),
           ("units", "eps_a_pct [%] compression positive | sr = sigma1'/sigma3' [-] | eps_v_pct [%] dilation negative")]
    rows_out = [(float(a), float(b), None) for a, b in ra] + [(float(a), None, float(b)) for a, b in ev]
    rows_out.sort(key=lambda r: (r[0], r[1] is None))
    os.makedirs(os.path.join(HERE, "tatsuoka1986"), exist_ok=True)
    with open(os.path.join(HERE, "tatsuoka1986", f"{tid}.csv"), "w", newline="", encoding="utf-8") as f:
        for k, v in hdr:
            f.write(f"# {k}: {v}\n")
        cw = csv.writer(f, lineterminator="\n")
        cw.writerow(["eps_a_pct", "sr", "eps_v_pct"])
        for a, b, c in rows_out:
            cw.writerow([fmt(a), fmt(b), fmt(c)])
    ip = int(np.argmax(ra[:, 1]))
    Rp = ra[ip, 1]
    phi = math.degrees(math.asin((Rp - 1) / (Rp + 1)))
    dil = None
    if len(ev):
        g = np.linspace(0, min(ev[-1, 0], 10), 401)
        evi = np.interp(g, ev[:, 0], ev[:, 1])
        k = 20
        dil = float(np.max(-(evi[k:] - evi[:-k]) / (g[k:] - g[:-k])))
    meta.append([tid, key, t["fig"], t["page"], t["s3"], round(t["s3"] / KGF, 4), t["e"],
                 round((E_MAX - t["e"]) / (E_MAX - E_MIN), 3), 90, len(ra), len(ev), fmt(Rp, 4), fmt(float(ra[ip, 0]), 3),
                 fmt(phi, 4), fmt(dil, 3), t["marker"]])
w(os.path.join(HERE, "tatsuoka1986", "tests_meta.csv"),
  ["Tatsuoka et al. 1986 test metadata. D_r_assumed = (0.977 - e_0.05)/(0.977 - 0.597): e_max/e_min are NOT in the paper (see README).",
   "R_peak, eps_a_at_peak, phi_peak and dilatancy_max are DERIVED here from the digitised curve (phi = asin((R-1)/(R+1)); dilatancy_max = max(-d eps_v / d eps_a) on a 0.5 % base, compression positive)."],
  ["test_id", "TIMs_key", "figure", "journal_page", "sigma3_kPa", "sigma3_kgf_cm2", "e_0.05", "D_r_assumed", "delta_deg",
   "n_ratio_pts", "n_epsv_pts", "R_peak", "eps_a_at_peak_pct", "phi_peak_deg_derived", "dilatancy_max_derived", "digitisation"], meta)

# FIG9 / FIG22 point series
rows9, rows22 = [], []
for s, pts in D["FIG9"].items():
    for e, ph in pts:
        rows9.append([float(s), round(float(s) * KGF, 2), e, ph])
for s, pts in D["FIG22"].items():
    for e, ep in pts:
        rows22.append([float(s), round(float(s) * KGF, 2), e, ep])
w(os.path.join(HERE, "tatsuoka1986", "fig9_phi_peak_vs_e.csv"),
  ["Tatsuoka et al. 1986 Fig. 9 (journal p. 75): phi_peak (method C-2-T, deg) vs e_0.05, delta = 90 deg, isotropically consolidated, plane strain.",
   "Digitised by TIMs (tests_data.py FIG9), copied unedited. Visual check against the PDF (this pack): all listed points agree to ~0.3 deg / 0.003 in e.",
   "KNOWN GAP: the sigma_c' = 0.1 kgf/cm2 series has only 2 points; the figure shows further 'x' markers (near e 0.746 and 0.805, phi ~44 and ~41 deg) that were not digitised."],
  ["sigma_c_kgf_cm2", "sigma_c_kPa", "e_0.05", "phi_peak_deg"], rows9)
w(os.path.join(HERE, "tatsuoka1986", "fig22_eps_peak_vs_e.csv"),
  ["Tatsuoka et al. 1986 Fig. 22 (journal p. 82): axial strain (external, boundary displacements) at (sigma1'/sigma3')max vs e_0.05, delta = 90 deg, plane strain.",
   "Digitised by TIMs (tests_data.py FIG22), copied unedited. Visual check against the PDF: points agree to ~0.1 % strain / 0.003 in e."],
  ["sigma_c_kgf_cm2", "sigma_c_kPa", "e_0.05", "eps_a_at_peak_pct"], rows22)

# joined point tests in the harness point-test schema: Fig. 9 rows, eps_peak from the Fig. 22 point at the same
# sigma_c' with |delta e| <= 0.0035 (global nearest-first, each point used once). Fig. 22 points without a Fig. 9 partner are not listed here.
jr = []
TOL = 0.0035
for s_, pts in D["FIG9"].items():
    cand = sorted(((abs(e9 - e22), i, j) for i, (e9, _) in enumerate(pts) for j, (e22, _) in enumerate(D["FIG22"][s_])))
    used9, used22, pair = set(), set(), {}
    for dist, i, j in cand:  # global greedy: smallest |delta e| first, each point used once
        if dist <= TOL and i not in used9 and j not in used22:
            used9.add(i)
            used22.add(j)
            pair[i] = D["FIG22"][s_][j][1]
    for i, (e, ph) in enumerate(pts):
        jr.append([round(float(s_) * KGF, 2), e, ph, fmt(pair.get(i)), "PS"])
w(os.path.join(HERE, "tatsuoka1986", "point_tests.csv"),
  ["Tatsuoka et al. 1986 point tests, delta = 90 deg, drained plane strain, e = e_0.05: phi_peak from Fig. 9 (p. 75), eps_peak from Fig. 22 (p. 82) where a point of the same sigma_c' lies within 0.0035 in e.",
   "Convenience join of fig9_phi_peak_vs_e.csv and fig22_eps_peak_vs_e.csv; unmatched eps_peak is left empty. Schema: harness/data.py load_points."],
  ["sigma3_kPa", "e", "phi_peak_deg", "eps_peak_pct", "kind"], jr)

# ---------------------------------------------------------------- Kimura 1985
KP = os.path.join(HERE, "kimura1985", "source", "kimura1985_fig9_digitized.csv")
lines = [l for l in open(KP, encoding="utf-8") if not l.startswith("#")]
rd = list(csv.DictReader(io.StringIO("".join(lines))))
by = {}
for r in rd:
    by.setdefault(r["test_id"], []).append(r)
km = []
GS = 2.64  # Tatsuoka 1986 p. 70 (Toyoura, specific gravity 2.64); Kimura gives none
for tid, rs in by.items():
    q = np.array([float(r["q_kPa"]) for r in rs])
    sB = np.array([float(r["s_over_B"]) for r in rs])
    S = np.array([float(r["S_model_mm"]) for r in rs])
    i = int(np.argmax(q))
    dr = float(tid[1:])
    e = E_MAX - dr / 100 * (E_MAX - E_MIN)
    km.append([tid, tid[0], dr, round(e, 4), round(GS * 9.81 / (1 + e), 2), 30, 30, 0.9, len(rs), fmt(float(q[i]), 5),
               fmt(float(sB[i]), 4), fmt(float(S[i]), 4), fmt(float(sB.max()), 4)])
w(os.path.join(HERE, "kimura1985", "tests_meta.csv"),
  ["Kimura, Kusakabe & Saitoh (1985) Geotechnique 35(1):33-45, Fig. 9 (journal p. 40): strip footing B = 30 mm, 30 g, Toyoura sand, pouring method, anisotropy series.",
   "Case V = load perpendicular to the bedding plane, H = parallel (p. 39, Fig. 8). D_r from the curve labels (Fig. 9). Prototype B = 30 mm x 30 = 0.90 m.",
   "e_assumed = 0.977 - D_r*(0.977-0.597) and gamma_derived = Gs*9.81/(1+e) with Gs = 2.64 are DERIVED/ASSUMED here: Kimura 1985 gives NO e_max, e_min, Gs or unit weight for the pouring series.",
   "q_peak, s_over_B_at_peak, S_at_peak and s_over_B_max are read from the digitised curve (source/kimura1985_fig9_digitized.csv)."],
  ["test_id", "case", "D_r_pct", "e_assumed", "gamma_derived_kN_m3", "B_model_mm", "g_level", "B_prototype_m", "n_points",
   "q_peak_kPa", "s_over_B_at_peak", "S_model_mm_at_peak", "s_over_B_max"], km)

# ---------------------------------------------------------------- DM04
T1 = [  # group, constant, symbol, value, role
    ("Elasticity", "G0", "G0", 125, "G = G0 p_at (2.97-e)^2/(1+e) (p/p_at)^0.5 (Table 2, p. 628)"),
    ("Elasticity", "nu", "nu", 0.05, "K = 2(1+nu)G/(3(1-2nu)) (Table 2)"),
    ("Critical state", "M", "M", 1.25, "triaxial compression CSL slope q/p"),
    ("Critical state", "c", "c", 0.712, "M_e/M_c, extension/compression ratio (Eq. 16)"),
    ("Critical state", "lambda_c", "lambda_c", 0.019, "CSL e_c = e0 - lambda_c (p_c/p_at)^xi (Table 2)"),
    ("Critical state", "e0", "e0", 0.934, "CSL intercept (e_c at p = 0), NOT the initial void ratio"),
    ("Critical state", "xi", "xi", 0.7, "CSL exponent"),
    ("Yield surface", "m", "m", 0.01, "yield-cone opening (Eq. 13)"),
    ("Plastic modulus", "h0", "h0", 7.05, "b0 = G0 h0 (1 - c_h e)(p/p_at)^-0.5 (Table 2)"),
    ("Plastic modulus", "ch", "c_h", 0.968, "idem"),
    ("Plastic modulus", "nb", "n^b", 1.1, "M^b = M exp(-n^b psi) (Eq. 9)"),
    ("Dilatancy", "A0", "A0", 0.704, "A_d = A0 (1 + <z:n>) (Eq. 27)"),
    ("Dilatancy", "nd", "n^d", 3.5, "M^d = M exp(n^d psi) (Eq. 10)"),
    ("Fabric-dilatancy", "zmax", "z_max", 4, "Eq. 26 (cyclic only)"),
    ("Fabric-dilatancy", "cz", "c_z", 600, "Eq. 26 (cyclic only)"),
]
w(os.path.join(HERE, "dm04", "toyoura_table1.csv"),
  ["Dafalias & Manzari (2004) J. Eng. Mech. 130(6):622-634, Table 1 'Model Constants' (PDF page 5 = journal p. 626); read visually from a 250-dpi render.",
   "Calibrated by DM04 to Verdugo & Ishihara (1996) triaxial data on Toyoura sand (p'=100-3000 kPa, D_r 18.5-63.7 %, e 0.735-0.907; text p. 632).",
   "p_at is NOT tabulated in DM04; G0 is dimensionless (multiplies p_at). Every dimensionless constant here matches the cluster's TOYOURA set and gate0_toyoura_oracle.py docstring."],
  ["group", "constant", "symbol", "value", "role_equation"], T1)

PAT = 100.0
LAM, E0C, XI = 0.019, 0.934, 0.7
rows = []
for p in (1.0, 4.9, 9.8, 49.0, 98.0, 100.0, 500.0, 1000.0, 3000.0):
    ec = E0C - LAM * (p / PAT) ** XI
    rows.append([p, round(ec, 4)] + [round(e - ec, 4) for e in (0.6517, 0.700, 0.714, 0.716, 0.755)])
w(os.path.join(HERE, "dm04", "toyoura_csl_derived.csv"),
  ["DERIVED arithmetic from DM04 Table 1: e_c(p) = 0.934 - 0.019 (p/p_at)^0.7 with p_at = 100 kPa (ASSUMED: not tabulated in DM04; the WP-134/150 oracle and the OpenSees example use 100, TIMs uses 101).",
   "psi = e - e_c at the listed e (state parameter). e = 0.6517 is Kimura's V85.6 under the 0.977/0.597 assumption; 0.700/0.714/0.716/0.755 are Tatsuoka's e_0.05 (void ratio at the test sigma_3' is lower and not reported)."],
  ["p_kPa", "e_c", "psi_at_e0.6517", "psi_at_e0.700", "psi_at_e0.714", "psi_at_e0.716", "psi_at_e0.755"], rows)
phi_c = math.degrees(math.asin(3 * 1.25 / (6 + 1.25)))
Me = 0.712 * 1.25
phi_e = math.degrees(math.asin(3 * Me / (6 - Me)))
print(f"DM04 phi_c = {phi_c:.3f} deg (M=1.25), phi_e = {phi_e:.3f} deg (M_e = c M = {Me:.3f})")

# ---------------------------------------------------------------- Ottawa F65
EMIN_OT, EMAX_OT = 0.4915, 0.7389  # averages of 9 trials, Characterization slide 7
zf = zipfile.ZipFile(ZIP)
names = sorted(n for n in zf.namelist() if n.startswith("Monotonic Triaxial Experiments/") and n.endswith(".txt"))
assert len(names) == 16, len(names)
cols_in = ["index", "time_min", "vertical_strain_pct", "volumetric_strain_pct", "corrected_area_mm2", "deviator_load_N",
           "deviator_stress_kPa", "pore_pressure_kPa", "horizontal_stress_kPa", "vertical_stress_kPa"]
om = []
for n in names:
    m = re.search(r"/(\d+)kPa_Drained_(Comp|Extension)_eo_0_(\d+)\.txt$", n)
    s_nom, mode, e3 = int(m.group(1)), ("TC" if m.group(2) == "Comp" else "TE"), m.group(3)
    e0 = float("0." + e3)
    txt = zf.read(n).decode("ascii", "replace").replace("\r", "")
    data = []
    for ln in txt.split("\n"):
        p = ln.split()
        if len(p) == 10 and re.match(r"^\d+$", p[0]):
            data.append([float(x) for x in p])
    a = np.array(data)
    assert a.shape[1] == 10 and len(a) > 50
    tid = f"ottawa_f65_{mode}_e{e3}_s{s_nom}"
    w(os.path.join(HERE, "ottawa_f65", f"{tid}.csv"),
      [f"Ottawa F65 sand, LEAP-2015, GWU (Vasko, El Ghoraiby, Manzari), consolidated drained triaxial {'compression' if mode == 'TC' else 'extension'}, e_0 = {e0}, nominal sigma_c' = {s_nom} kPa",
       f"source: download.zip :: {n}   (values copied verbatim; only the 3 header rows were replaced by this one)",
       "units: time min | strains % | area mm^2 | load N | stresses kPa (TOTAL stresses; effective = total - pore_pressure)",
       "sign: compression positive. vertical_strain < 0 and deviator < 0 in extension; volumetric_strain < 0 = dilation (matches the volumetric plots of Monotonic Test Results.pdf)"],
      cols_in, [[fmt(x, 8) if i else int(x) for i, x in enumerate(r)] for r in a])
    v_eff = a[:, 9] - a[:, 7]
    h_eff = a[:, 8] - a[:, 7]
    s1 = np.maximum(v_eff, h_eff)
    s3 = np.minimum(v_eff, h_eff)
    ph = np.degrees(np.arcsin((s1 - s3) / (s1 + s3)))
    ip = int(np.argmax(ph))
    q = a[:, 6]
    iq = int(np.argmax(np.abs(q)))
    p0 = (a[0, 9] + 2 * a[0, 8]) / 3 - a[0, 7]
    om.append([tid, mode, e0, round((EMAX_OT - e0) / (EMAX_OT - EMIN_OT), 3), s_nom, round(float(h_eff[0]), 2), round(float(p0), 2),
               len(a), fmt(float(a[-1, 2]), 4), fmt(float(q[iq]), 5), fmt(float(a[iq, 2]), 4), fmt(float(ph[ip]), 4),
               fmt(float(a[ip, 2]), 4), fmt(float(a[-1, 3]), 4), n.replace("Monotonic Triaxial Experiments/", "")])
w(os.path.join(HERE, "ottawa_f65", "tests_meta.csv"),
  ["Ottawa F65, 16 monotonic drained triaxial tests (7 loose e0 0.724: 6 TC + 1 TE; 3 dense e0 0.584 TC; 6 'test density' e0 0.604: 3 TC + 3 TE).",
   "D_r_derived = (0.7389 - e0)/(0.7389 - 0.4915): averages of 9 e_max/e_min trials, Characterization Tests.pdf slide 7 (std 0.0247 / 0.0183 -> D_r at e0 0.724 is within the scatter of e_max).",
   "sigma3_eff_start = (horizontal - pore) of row 1; p_eff_start from row 1. q_peak_signed = deviator stress at max |q|. phi_peak_deg_derived = max asin((s1'-s3')/(s1'+s3')) from effective stresses (s1' = larger of vertical/horizontal).",
   "eps_a_at_phi_peak and eps_v_end are read from the file (compression positive)."],
  ["test_id", "mode", "e0", "D_r_derived", "sigma3_nominal_kPa", "sigma3_eff_start_kPa", "p_eff_start_kPa", "n_rows", "eps_a_end_pct",
   "q_peak_signed_kPa", "eps_a_at_q_peak_pct", "phi_peak_deg_derived", "eps_a_at_phi_peak_pct", "eps_v_end_pct", "source_file"], om)
print("done:", len(names), "ottawa tests;", len(meta), "tatsuoka tests;", len(km), "kimura curves")
