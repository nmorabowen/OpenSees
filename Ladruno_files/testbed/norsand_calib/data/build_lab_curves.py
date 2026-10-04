"""WP-144 P3 data pack: lab-curve converters for the Wang Toyoura database and the Ottawa F65 tests.

    python Ladruno_files/testbed/norsand_calib/data/build_lab_curves.py <Wang_sand_triaxial_database_rev2.zip>

Reads    the Wang database zip (read only; never modified): integrated_dataset/33_Toyoura_sand_HKU.csv, 34_Toyoura_sand_Tokyo.csv
         data/ottawa_f65/ottawa_f65_*.csv   (the verbatim LEAP-2015 files already in the pack)
Writes   data/wang_toyoura/wang{33_HKU|34_Tokyo}_TMD<n>.csv   one harness lab-curve file per test (harness/data.py schema)
         data/wang_toyoura/tests_meta.csv
         data/ottawa_f65/curves/<test_id>.csv                  harness lab-curve files, sr from EFFECTIVE stresses
         data/ottawa_f65/curves/curves_meta.csv
Nothing here fits or tunes anything. Derived columns are labelled "derived".

Conventions (harness/data.py):
  eps_a_pct   as in the source file; the loader multiplies it by eps_a_sign. The harness peak logic (LabCurve.peak, x_end = eps_peak + post_peak)
              needs the LOADED axial strain to be POSITIVE, so TE files (Ottawa vertical_strain < 0 in extension) carry eps_a_sign = -1.
  sr          sigma1'/sigma3' = major/minor effective principal stress (TE: horizontal/vertical, > 1).
  eps_v_pct   volumetric strain %, compression positive, dilation NEGATIVE ('dilation_negative'); both Wang and Ottawa are already so.
"""
import csv
import math
import os
import sys
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ZIP = sys.argv[1]


def fmt(x, n=7):
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{n}g}"


def write_curve(path, meta, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        for k, v in meta.items():
            f.write(f"# {k}: {v}\n")
        cw = csv.writer(f, lineterminator="\n")
        cw.writerow(["eps_a_pct", "sr", "eps_v_pct"])
        for r in rows:
            cw.writerow([fmt(float(x)) for x in r])


def phi_of(sr):
    return math.degrees(math.asin((sr - 1.0) / (sr + 1.0)))


# ------------------------------------------------------------------------------------------- Wang (HKU, Tokyo)
def read_wang(text):
    sec, hdr, out = None, None, {}
    for line in text.splitlines():
        if line.startswith("# ====="):
            sec = line.split("]")[1].strip(" =")
            out[sec], hdr = [], None
        elif line.startswith("#") or not line.strip() or sec is None:
            continue
        else:
            r = next(csv.reader([line]))
            if hdr is None:
                hdr = r
            else:
                out[sec].append(dict(zip(hdr, r)))
    return out


WANG = [("33_Toyoura_sand_HKU", "wang33_HKU", "Chen & Yang (2025) Eng. Geol. 345:107863 (University of Hong Kong)"),
        ("34_Toyoura_sand_Tokyo", "wang34_Tokyo", "Verdugo & Ishihara (1996) Soils and Foundations 36(2):81-91 (University of Tokyo)")]
os.makedirs(os.path.join(HERE, "wang_toyoura"), exist_ok=True)
zf = zipfile.ZipFile(ZIP)
meta_rows = []
for stem, tag, src in WANG:
    d = read_wang(zf.read(f"integrated_dataset/{stem}.csv").decode("utf-8"))
    idx = {r["symbol"]: r["value"] for r in d["INDEX_PROPERTIES"]}
    for t in d["TEST_PROGRAMME"]:
        rows = [r for r in d["TEST_DATA"] if r["test_id"] == t["test_id"]]
        ea = np.array([float(r["eps_a[%]"]) for r in rows])
        ev = np.array([float(r["eps_v[%]"]) if r["eps_v[%]"] else float("nan") for r in rows])
        sr_ = np.array([float(r["sigma_r[kPa]"]) for r in rows])
        sv = np.array([float(r["sigma_v[kPa]"]) for r in rows])
        R = sv / sr_
        ip = int(np.argmax(R))
        tid = f"{tag}_{t['test_id']}"
        s3 = float(t["sigma_r0[kPa]"])
        e0 = float(t["e_0[-]"])
        ev0 = float(ev[0]) if len(ev) else float("nan")
        peak_at_end = ip >= len(R) - 3
        meta = {
            "id": tid, "kind": "TC", "sigma3_kPa": s3, "e0": e0,
            "e0_definition": "e_0 of the Wang database TEST_PROGRAMME = void ratio after isotropic consolidation, at the start of shearing (as reported by the source)",
            "source": f"Wang sand triaxial database rev2 (4TU, doi 10.4121/086847a6-ba39-4d66-973b-6b93028c7ad8, CC-BY-4.0), file {stem}.csv, {t['test_id']}; original data: {src}",
            "eps_v_convention": "dilation_negative", "eps_a_sign": 1,
            "material": f"Toyoura sand ({'HKU' if 'HKU' in stem else 'Tokyo'}), Gs {idx['Gs']}, e_max {idx['e_max']}, e_min {idx['e_min']} (the database's own values for this file)",
            "test": "drained isotropically consolidated triaxial compression; sigma_r held; sr = sigma_v'/sigma_r' (effective stresses as given)",
            "Dr_0_pct_database": t["Dr_0[%]"],
            "digitisation": "Wang compiled these from published FIGURES (digitised; strains rounded to 4 decimals, stresses to 0.01 kPa); tolerance not stated by the database",
            "eps_v_at_eps_a0_pct": fmt(ev0, 4),
            "peak_at_end_of_test": str(peak_at_end).lower(),
            "units": "eps_a_pct [%] compression positive | sr [-] | eps_v_pct [%] dilation negative",
        }
        if t["notes"]:
            meta["database_note"] = t["notes"]
        write_curve(os.path.join(HERE, "wang_toyoura", f"{tid}.csv"), meta, list(zip(ea, R, ev)))
        # dilatancy: max over a 1 % base of -d eps_v / d eps_a
        g = np.arange(0, ea[-1], 0.1)
        evi = np.interp(g, ea, ev)
        k = 10
        dil = float(np.max(-(evi[k:] - evi[:-k]) / (g[k:] - g[:-k]))) if len(g) > k else float("nan")
        meta_rows.append([tid, stem, t["test_id"], s3, e0, t["Dr_0[%]"], len(rows), fmt(float(ea[0]), 3), fmt(float(ea[-1]), 4),
                          fmt(float(R[ip]), 4), fmt(float(ea[ip]), 4), fmt(phi_of(float(R[ip])), 4), str(peak_at_end).lower(),
                          fmt(float(ev[-1]), 4), fmt(ev0, 4), fmt(dil, 3), t["notes"]])
with open(os.path.join(HERE, "wang_toyoura", "tests_meta.csv"), "w", newline="", encoding="utf-8") as f:
    f.write("# Wang database Toyoura tests (33 HKU, 34 Tokyo) in the harness lab-curve schema. R_peak, eps_a_at_peak, phi_peak_deg_derived, dilatancy_max_derived are DERIVED here from the curve\n")
    f.write("# (phi = asin((R-1)/(R+1)); dilatancy_max = max(-d eps_v/d eps_a) on a 1 % base). peak_at_end = max R within the last 3 points (no interior peak: contractive/loose tests).\n")
    cw = csv.writer(f, lineterminator="\n")
    cw.writerow(["test_id", "database_file", "database_test_id", "sigma3_kPa", "e0", "Dr0_pct_database", "n_points", "eps_a_first_pct",
                 "eps_a_end_pct", "R_peak", "eps_a_at_peak_pct", "phi_peak_deg_derived", "peak_at_end", "eps_v_end_pct",
                 "eps_v_at_eps_a0_pct", "dilatancy_max_derived", "database_note"])
    cw.writerows(meta_rows)

# ------------------------------------------------------------------------------------------- Ottawa F65
OD = os.path.join(HERE, "ottawa_f65")
os.makedirs(os.path.join(OD, "curves"), exist_ok=True)
with open(os.path.join(OD, "tests_meta.csv"), encoding="utf-8") as f:
    om = list(csv.DictReader(l for l in f if not l.startswith("#")))
orows = []
for m in om:
    tid = m["test_id"]
    with open(os.path.join(OD, tid + ".csv"), encoding="utf-8") as f:
        rd = list(csv.DictReader(l for l in f if not l.startswith("#")))
    ea = np.array([float(r["vertical_strain_pct"]) for r in rd])
    ev = np.array([float(r["volumetric_strain_pct"]) for r in rd])
    u = np.array([float(r["pore_pressure_kPa"]) for r in rd])
    sh = np.array([float(r["horizontal_stress_kPa"]) for r in rd]) - u   # effective radial stress
    sv = np.array([float(r["vertical_stress_kPa"]) for r in rd]) - u     # effective axial stress
    sr = np.maximum(sv, sh) / np.minimum(sv, sh)                         # sigma1'/sigma3' (TC: sv/sh ; TE: sh/sv)
    te = m["mode"] == "TE"
    ip = int(np.argmax(sr))
    meta = {
        "id": tid, "kind": m["mode"], "sigma3_kPa": float(m["sigma3_nominal_kPa"]), "e0": float(m["e0"]),
        "e0_definition": "e_0 label of the LEAP-2015 test (rho_d from the file/folder name with Gs 2.65; see data/README.md section 6 item 8b); void ratio at the start of shearing",
        "source": "George Washington University LEAP-2015 monotonic triaxial database (Vasko, El Ghoraiby, Manzari, 2014), file " + m["source_file"],
        "eps_v_convention": "dilation_negative", "eps_a_sign": -1 if te else 1,
        "eps_a_sign_note": ("the file's vertical_strain is NEGATIVE in extension; -1 makes the loaded eps_a positive, which harness LabCurve.peak / objective x_end need"
                            if te else "compression positive in the file"),
        "material": "Ottawa F65 sand, saturated, drained, isotropically consolidated",
        "test": "drained triaxial " + ("EXTENSION (axial unloading, radial stress held): sr = sigma_h'/sigma_v' = sigma1'/sigma3'" if te else "compression: sr = sigma_v'/sigma_h'"),
        "stress_method": "sr from EFFECTIVE stresses = (total - pore_pressure_kPa) per row; sigma3_kPa is the nominal consolidation stress (effective start " + m["sigma3_eff_start_kPa"] + " kPa)",
        "Dr_derived": m["D_r_derived"],
        "peak_at_end_of_test": str(ip >= len(sr) - 3).lower(),
        "units": "eps_a_pct [%] as in the file (TE negative; times eps_a_sign) | sr [-] | eps_v_pct [%] dilation negative",
    }
    if tid == "ottawa_f65_TE_e604_s100":
        meta["caveat"] = ("anomalous: phi_peak 44.1 deg vs 35.0/36.4 at 200/300 kPa, same density; q only ~-80 kPa (low end of load cell); "
                          "82 rows, stops at -8.1 %. Low weight or exclude (data/README.md 6.8c)")
    write_curve(os.path.join(OD, "curves", tid + ".csv"), meta, list(zip(ea, sr, ev)))
    orows.append([tid, m["mode"], m["e0"], m["sigma3_nominal_kPa"], len(rd), fmt(float(ea[0]), 3), fmt(float(ea[-1]), 4), fmt(float(sr[ip]), 4),
                  fmt(float(ea[ip]), 4), fmt(phi_of(float(sr[ip])), 4), str(ip >= len(sr) - 3).lower(), -1 if te else 1, fmt(float(ev[-1]), 4)])
with open(os.path.join(OD, "curves", "curves_meta.csv"), "w", newline="", encoding="utf-8") as f:
    f.write("# Ottawa F65 harness lab curves (converter build_lab_curves.py): sr from effective stresses; eps_a_sign -1 for TE. R_peak, phi derived from the converted sr.\n")
    cw = csv.writer(f, lineterminator="\n")
    cw.writerow(["test_id", "kind", "e0", "sigma3_nominal_kPa", "n_points", "eps_a_first_pct", "eps_a_end_pct", "R_peak", "eps_a_at_peak_pct",
                 "phi_peak_deg_derived", "peak_at_end", "eps_a_sign", "eps_v_end_pct"])
    cw.writerows(orows)
print("wang:", len(meta_rows), "tests; ottawa:", len(orows), "tests")
