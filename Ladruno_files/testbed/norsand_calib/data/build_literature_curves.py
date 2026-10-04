"""WP-144 P3 data pack: literature stress-strain curves (Fukushima & Tatsuoka 1984 Fig. 5, Lam & Tatsuoka 1988 Figs. 5 and 9) -> harness lab-curve files.

    python Ladruno_files/testbed/norsand_calib/data/build_literature_curves.py

Reads    fukushima1984/digitised_points_fig5.json, fukushima1984/table1_tests.csv, lam_tatsuoka1988/digitised_points_fig5_fig9.json
Writes   fukushima1984/curves/<id>.csv + curves_meta.csv ; lam_tatsuoka1988/curves/<id>.csv + curves_meta.csv   (harness/data.py lab-curve schema)
The digitised points (JSON) are the output of the digitising pipeline described in data/README.md section 8 (tools kept in data/_digitise/).
Nothing here fits anything.

Strain conversions (compression positive; eps_v dilation NEGATIVE):
  Fukushima & Tatsuoka: the abscissa is the axial strain eps_a itself.
  Lam & Tatsuoka TC  : abscissa D = eps1 - eps3 = (3 eps1 - eps_v)/2          ->  eps_a = eps1 = (2 D + eps_v)/3
  Lam & Tatsuoka PSC : abscissa D = eps1 - eps3, eps2 = 0, eps_v = eps1+eps3   ->  eps_a = eps1 = (D + eps_v)/2
  eps_v at the stress-ratio abscissae is interpolated from the same test's eps_v trace (origin (0,0) added); stress-ratio rows with D beyond
  the eps_v trace by more than 3 % are dropped; within 3 % eps_v is extrapolated linearly from the last segment.
"""
import csv
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KG = 98.0


def phi_of(sr):
    return math.degrees(math.asin((sr - 1.0) / (sr + 1.0)))


def write_curve(path, meta, rows, extra_cols=()):
    with open(path, "w", newline="", encoding="utf-8") as f:
        for k, v in meta.items():
            f.write(f"# {k}: {v}\n")
        cw = csv.writer(f, lineterminator="\n")
        cw.writerow(["eps_a_pct", "sr", "eps_v_pct", *extra_cols])
        for r in rows:
            cw.writerow(["" if x is None else (f"{x:.5g}" if isinstance(x, float) else x) for x in r])


def interp_ev(Dq, ev_pts, kmax=3.0):
    pts = sorted([(0.0, 0.0)] + [tuple(p) for p in ev_pts])
    x = np.array([p[0] for p in pts])
    y = np.array([p[1] for p in pts])
    out = []
    for D in Dq:
        if D <= x[-1]:
            out.append(float(np.interp(D, x, y)))
        elif D <= x[-1] + kmax and len(x) >= 3:
            s = (y[-1] - y[-3]) / (x[-1] - x[-3])
            out.append(float(y[-1] + s * (D - x[-1])))
        else:
            out.append(None)
    return out


# ----------------------------------------------------------------------------------------- Fukushima & Tatsuoka 1984 Fig. 5
FT = json.load(open(os.path.join(HERE, "fukushima1984", "digitised_points_fig5.json")))
with open(os.path.join(HERE, "fukushima1984", "table1_tests.csv"), encoding="utf-8") as f:
    T1 = {int(r["test"]): r for r in csv.DictReader(line for line in f if not line.startswith("#"))}
# panel, sigma_c' [kgf/cm2] -> (Table 1 test number, e0.3 from the figure legend)
FTMAP = {
    "5a": {"0.1": (28, 0.650), "0.2": (34, 0.658), "0.5": (45, 0.658), "1.0": (53, 0.671), "2.0": (66, 0.677), "4.0": (73, 0.679)},
    "5b": {"0.02": (9, 0.679), "0.05": (19, 0.681), "0.1": (29, 0.687)},
    "5d": {"0.02": (11, 0.853), "0.05": (21, 0.898), "0.1": (33, 0.866)},
}
PAGE = {
    "5a": "Fig. 5(a), journal p. 36 (PDF page 7), dense, sigma_c' 0.1-4 kgf/cm2",
    "5b": "Fig. 5(b), journal p. 36 (PDF page 7), dense, sigma_c' 0.02-0.1 kgf/cm2",
    "5d": "Fig. 5(d), journal p. 36 (PDF page 7), loose, sigma_c' 0.02-0.1 kgf/cm2",
}
outdir = os.path.join(HERE, "fukushima1984", "curves")
os.makedirs(outdir, exist_ok=True)
meta_rows = []
for panel, mp in FTMAP.items():
    for sc, (tno, e03) in mp.items():
        key = f"{panel}_{sc}"
        d = FT[key]
        t1 = T1[tno]
        sr = [tuple(p) for p in d["sr"]]
        ev = [tuple(p) for p in d["ev"]]
        srmax = max(y for x, y in sr)
        xpk = [x for x, y in sr if y == srmax][0]
        sr_t1 = float(t1["sr_derived"])
        phi_t1 = float(t1["phi_deg_uncorrected"])
        sigc = float(sc) * KG
        tid = f"ft84_fig{panel}_sc{sc.replace('.', 'p')}"
        s3f = t1["sigma3_f_kPa"]
        meta = {
            "id": tid,
            "kind": "TC",
            "sigma3_kPa": round(sigc, 2),
            "e0": e03,
            "e0_definition": "e_0.3 = void ratio at sigma_c' = 0.3 kgf/cm2 (29.4 kPa), the paper's density label (figure legend); for sigma_c' != 0.3 the void ratio at the shear stress is NOT reported (p. 35, Fig. 4)",
            "source": f"Fukushima & Tatsuoka (1984) Soils and Foundations 24(4):30-48, {PAGE[panel]}; Table 1 test {tno} (p. 32)",
            "eps_v_convention": "dilation_negative",
            "eps_a_sign": 1,
            "phi_peak_deg": phi_t1,
            "eps_peak_pct": round(xpk, 2),
            "phi_peak_note": f"Table 1 phi (mid-height, UNCORRECTED for membrane forces, t0 = {t1['t0_mm']} mm); the digitised maximum sr {srmax:.3f} vs {sr_t1:.3f} from Table 1 (difference {srmax - sr_t1:+.3f}) is the calibration check of the digitisation",
            "sigma3_f_kPa": s3f if s3f != "" else "not measured (sigma_c' = 4.0 kgf/cm2)",
            "material": "Toyoura sand (D50 0.16 mm, Uc 1.46, Gs 2.64, e_max 0.977, e_min 0.605), saturated, air-pluviated, 7 cm diameter x 15 cm, lubricated ends (p. 33)",
            "test": "drained isotropically consolidated triaxial compression, axial strain rate 0.25 %/min; sr = sigma1'/sigma3' at the MID-HEIGHT of the sample, UNCORRECTED for membrane forces (the figure title says so); eps_a = axial strain",
            "caveat_membrane": f"membrane t0 = {t1['t0_mm']} mm, latex, E_m 15.2 kgf/cm2: at sigma_c' <= 0.1 kgf/cm2 the uncorrected stress ratios are biased (paper Figs. 7, 12-17); at the lowest stresses sigma3' drifted during the test: (sigma3')_f = {s3f} kPa vs sigma_c' = {sigc:.2f} kPa",
            "uncertainty": "sr +-0.05 (1 px = 0.01; calibration cross-checked against Table 1 phi, see phi_peak_note), +-0.1 where markers of neighbouring curves merge; eps_a +-0.1 %; eps_v +-0.15 % (scan skew handled by a bilinear map between the frame lines)",
            "digitisation": "marker-chain / column-cluster tracker on the native 1-bit scan (data/_digitise); rows fill ONE of sr / eps_v (the two traces are digitised at different eps_a); the first row (0,1,0) is the assumed isotropic start, not digitised",
            "units": "eps_a_pct [%] compression positive | sr [-] | eps_v_pct [%] dilation negative",
        }
        if d.get("note"):
            meta["trace_note"] = d["note"]
        rows = [(0.0, 1.0, 0.0)] + [(x, y, None) for x, y in sr] + [(x, None, y) for x, y in ev]
        rows.sort(key=lambda r: (r[0], r[1] is None))
        write_curve(os.path.join(outdir, tid + ".csv"), meta, rows)
        meta_rows.append([
            tid, panel, sc, round(sigc, 2), e03, tno, t1["t0_mm"], len(sr), len(ev),
            f"{min(x for x, y in sr):.2f}", f"{max(x for x, y in sr):.2f}",
            f"{min(x for x, y in ev):.2f}" if ev else "", f"{max(x for x, y in ev):.2f}" if ev else "",
            f"{srmax:.3f}", f"{xpk:.2f}", f"{sr_t1:.3f}", f"{srmax - sr_t1:+.3f}", f"{phi_of(srmax):.2f}", phi_t1,
            f"{min(y for x, y in ev):.2f}" if ev else "",
        ])
with open(os.path.join(outdir, "curves_meta.csv"), "w", newline="", encoding="utf-8") as f:
    f.write("# Fukushima & Tatsuoka (1984) Fig. 5 digitised curves.  sr_table1 = (1+sin phi)/(1-sin phi) from Table 1 (uncorrected phi); dsr = digitised max - Table 1 (calibration check).  phi_dig derived here.\n")
    cw = csv.writer(f, lineterminator="\n")
    cw.writerow(["id", "panel", "sigma_c_kgf_cm2", "sigma_c_kPa", "e0p3", "table1_test", "t0_mm", "n_sr", "n_ev", "ea_sr_first", "ea_sr_last",
                 "ea_ev_first", "ea_ev_last", "sr_max_digitised", "ea_at_sr_max", "sr_table1", "dsr", "phi_dig_deg", "phi_table1_deg", "ev_min_digitised"])
    cw.writerows(meta_rows)
print("FT curves:", len(meta_rows))

# ----------------------------------------------------------------------------------------- Lam & Tatsuoka 1988 Figs. 5 and 9
LT = json.load(open(os.path.join(HERE, "lam_tatsuoka1988", "digitised_points_fig5_fig9.json")))
outdir = os.path.join(HERE, "lam_tatsuoka1988", "curves")
os.makedirs(outdir, exist_ok=True)
rows_meta = []


def lt_meta(tid, kind, e03, source, test, extra):
    m = {
        "id": tid, "kind": kind, "sigma3_kPa": 98.0, "e0": e03,
        "e0_definition": "e_0.3 = void ratio at sigma3' = 0.3 kgf/cm2 (29.4 kPa), the paper's density label (legend); the specimen was then consolidated to 1.0 kgf/cm2 = 98 kPa and the void ratio at the start of shear is NOT reported (lower than e_0.3 by an unknown amount)",
        "source": source, "eps_v_convention": "dilation_negative", "eps_a_sign": 1,
        "material": "Toyoura sand (D50 0.16 mm, Uc 1.46, Gs 2.64, e_max 0.977, e_min 0.605), saturated, air-pluviated, prismatic specimens, lubricated rigid planes + flexible planes (mixed boundaries, Type 4 lubrication)",
        "test": test,
    }
    m.update(extra)
    return m


# ---- TC (Fig. 5)
TCE = {"a": {"0": 0.671, "30": 0.637, "60": 0.648, "90": 0.675}, "b": {"0": 0.657, "30": 0.649, "60": 0.642, "90": 0.671}}
HW = {"a": 1.0, "b": 0.25}
for pnl in ("a", "b"):
    for om in ("0", "30", "60", "90"):
        evk = f"5{pnl}_ev_{om}"
        srk = f"5{pnl}_sr_{om}" if f"5{pnl}_sr_{om}" in LT else f"5{pnl}_sr_6090"
        merged = srk.endswith("6090")
        ev = [(p[0], -2.0 * p[1]) for p in LT[evk]]      # stored ordinate is sr-equivalent y; eps_v = -2 y (right axis)
        sr = [tuple(p) for p in LT[srk]]
        Dq = [x for x, y in sr]
        evq = interp_ev(Dq, ev)
        rowsr = [((2 * D + e) / 3.0, y, None) for (D, y), e in zip(sr, evq) if e is not None]
        rowev = [((2 * D + e) / 3.0, None, e) for D, e in ev]
        rows = [(0.0, 1.0, 0.0)] + rowsr + rowev
        rows.sort(key=lambda r: (r[0], r[1] is None))
        srmax = max(r[1] for r in rowsr)
        xpk = [r[0] for r in rowsr if r[1] == srmax][0]
        tid = f"lt88_fig5{pnl}_w{int(om):02d}_hw{str(HW[pnl]).replace('.', 'p')}"
        extra = {
            "omega_deg": om, "xi_deg": "n/a (TC: sigma2 = sigma3)", "H_over_W": HW[pnl],
            "x_conversion": "paper abscissa D = eps1 - eps3 = (3 eps1 - eps_v)/2  ->  eps_a = (2 D + eps_v)/3 with eps_v(D) interpolated from this test's eps_v trace",
            "uncertainty": "sr +-0.04 (62 px per unit); D +-0.1 % (31 px per %), eps_v +-0.1 %, eps_a +-0.15 % after conversion; the initial rise (D < ~2 %) of sr is NOT resolved (the four curves overlap): the first row is the assumed isotropic start",
            "caveat_bedding_error": "strains not corrected for bedding error at the lubrication layers (text p. 94): eps_a overestimated; the relative error is larger for H/W = 0.25",
            "digitisation": "column-cluster family tracker on the native scan (data/_digitise/run_lt5.py); rows fill ONE of sr / eps_v",
            "units": "eps_a_pct [%] compression positive | sr [-] | eps_v_pct [%] dilation negative",
        }
        if pnl == "b":
            extra["caveat_H_over_W"] = "H/W = 0.25: strong end restraint (paper Figs. 4(b), 8(b)): phi larger than for H/W 1.0-2.0; use for the omega/anisotropy effect, not as an element test"
        if merged:
            extra["sr_note"] = "the sigma1'/sigma3' curves of omega = 60 and 90 coincide for D > ~5 % in this panel: ONE merged trace is given for both files"
        meta = lt_meta(
            tid, "TC", TCE[pnl][om],
            f"Lam & Tatsuoka (1988) Soils and Foundations 28(1):89-106, Fig. 5({pnl}), journal p. 94 (PDF page 6), TC-3, H/W = {HW[pnl]}, omega = {om} deg",
            "drained isotropically consolidated triaxial compression (TC-3), sigma3' = 98 kPa held, axial strain rate 0.25 %/min; sr = sigma1'/sigma3'", extra)
        meta["phi_peak_deg"] = round(phi_of(srmax), 2)
        meta["eps_peak_pct"] = round(xpk, 2)
        write_curve(os.path.join(outdir, tid + ".csv"), meta, rows)
        rows_meta.append([tid, "TC", TCE[pnl][om], 98.0, HW[pnl], om, len(rowsr), len(rowev), f"{srmax:.3f}", f"{phi_of(srmax):.2f}", f"{xpk:.2f}", f"{min(r[2] for r in rowev):.2f}"])

# ---- PSC (Fig. 9)
band = [(D, -2.0 * y) for D, y in LT["9_band_ev_y"]]          # eps_v of the overlapped triangle+circle band, D <= 8
PSC = {"tri": (1.9, 0.656), "ci": (1.0, 0.652)}
for key, (hw, e03) in PSC.items():
    sr = [tuple(p) for p in LT[f"9_sr_{key}"]]
    s2 = [tuple(p) for p in LT[f"9_s2_{key}"]]
    evm = [(x, -2.0 * y) for x, y in LT[f"9_ev_{key}"]]
    ev = sorted([p for p in band if p[0] < evm[0][0] - 0.2] + evm)
    Dq = [x for x, y in sr]
    evq = interp_ev(Dq, ev)
    rowsr = [((D + e) / 2.0, y, None, None) for (D, y), e in zip(sr, evq) if e is not None]
    s2e = interp_ev([x for x, y in s2], ev)
    rows2 = [((D + e) / 2.0, None, None, y) for (D, y), e in zip(s2, s2e) if e is not None]
    rowev = [((D + e) / 2.0, None, e, None) for D, e in ev]
    rows = [(0.0, 1.0, 0.0, 1.0)] + rowsr + rows2 + rowev
    rows.sort(key=lambda r: (r[0], r[1] is None))
    srmax = max(r[1] for r in rowsr)
    xpk = [r[0] for r in rowsr if r[1] == srmax][0]
    tid = f"lt88_fig9_psc_w00_hw{str(hw).replace('.', 'p')}"
    extra = {
        "omega_deg": 0, "xi_deg": 90, "H_over_W": hw, "b_at_failure_range": "0.20-0.34 (Fig. 8(a) lower lines, not transcribed)",
        "x_conversion": "paper abscissa D = eps1 - eps3; plane strain (eps2 = 0): eps_v = eps1 + eps3 -> eps_a = eps1 = (D + eps_v)/2; eps_v(D) for D < 8 % is the centre line of the overlapped triangle+circle eps_v band (both tests, +-0.15 %), for D >= 8 % the test's own markers",
        "uncertainty": "sr and sigma2'/sigma3' +-0.05 (marker centres; 64 px per unit); D +-0.15 % (32 px per %); eps_v +-0.15 % (band) / +-0.1 % (markers); eps_a +-0.2 % after conversion",
        "gaps": "only markers that are individually resolved are read: sigma1'/sigma3' from D ~1.4 % (initial rise not resolved); sigma2'/sigma3' only 4 markers (the triangle and circle curves overlap in a thick band for D = 6-11 %); eps_v band + markers; the square (H/W 0.5) and cross (H/W 0.25) series are NOT digitised (strong end restraint)",
        "extra_columns": "s2_over_s3 = sigma2'/sigma3' (rows with only this column filled)",
        "units": "eps_a_pct [%] compression positive | sr [-] = sigma1'/sigma3' | eps_v_pct [%] dilation negative | s2_over_s3 [-]",
    }
    meta = lt_meta(
        tid, "PS", e03,
        f"Lam & Tatsuoka (1988) Soils and Foundations 28(1):89-106, Fig. 9, journal p. 97 (PDF page 8), PSC-3, omega = 0, H/W = {hw}",
        "drained isotropically consolidated plane strain compression (PSC-3), sigma3' = 98 kPa held, eps2 = 0 (rigid lubricated planes), axial strain rate 0.25 %/min", extra)
    meta["phi_peak_deg"] = round(phi_of(srmax), 2)
    meta["eps_peak_pct"] = round(xpk, 2)
    write_curve(os.path.join(outdir, tid + ".csv"), meta, rows, extra_cols=("s2_over_s3",))
    rows_meta.append([tid, "PS", e03, 98.0, hw, "0", len(rowsr), len(rowev), f"{srmax:.3f}", f"{phi_of(srmax):.2f}", f"{xpk:.2f}", f"{min(r[2] for r in rowev):.2f}"])

with open(os.path.join(outdir, "curves_meta.csv"), "w", newline="", encoding="utf-8") as f:
    f.write("# Lam & Tatsuoka (1988) digitised curves (Figs. 5, 9).  phi and eps at peak are derived from the digitised sr (marker resolution: the true maximum may lie between markers).\n")
    cw = csv.writer(f, lineterminator="\n")
    cw.writerow(["id", "kind", "e0p3", "sigma3_kPa", "H_over_W", "omega_deg", "n_sr", "n_ev", "sr_max_digitised", "phi_dig_deg", "ea_at_sr_max_pct", "ev_min_pct"])
    cw.writerows(rows_meta)
print("LT curves:", len(rows_meta))
