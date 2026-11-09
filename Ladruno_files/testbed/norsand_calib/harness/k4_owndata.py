"""K4, own-data form (owner decision O8, 2026-11-08): each model against the data it was calibrated on.

The cross-model K4 (k4_crossmodel) compares NorSand with DM04 `uw_model`; the two were fitted to different data
(NorSand: the K3 set, mostly plane strain; DM04: Verdugo & Ishihara 1996 drained triaxial compression), so a
difference between them mixes model form with calibration. Here each model meets data, with the same metric and the
same tolerance as the cross-model gate (K3's scatter model, objective.Weights: |d sr_peak| <= 0.10 and
|d ev_peak| <= 0.10 %, ev at the peak point, compression positive):

                      K3 set (19 curves + 22 point tests)    Verdugo & Ishihara 1996 (6 drained TC curves)
  NorSand (K3 theta)  GATED                                   reported
  DM04 uw_model       reported                                GATED

The off-diagonal cells say how much of the cross-model gap is calibration data. Data: data/README.md (Tatsuoka
1986, Lam & Tatsuoka 1988 omega 0, Fukushima & Tatsuoka 1984 Fig. 5, exactly fit_tatsuoka --lit; V&I = the Wang
database file 34_Toyoura_sand_Tokyo, data/wang_toyoura/wang34_Tokyo_*). Point tests carry phi_peak only (sr gated,
ev not available); eps_peak is reported. A data curve without an eps_v branch is gated on sr alone.
p_a is TIMS_P_A (101 kPa) for both models, as in k4_crossmodel; DM04's own Toyoura runs used 100 (open item O9).

  python -m harness.k4_owndata --theta <K3 result JSON | JSON dict> --out out/k4_own [--workers 8]
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from . import data as DATA
from .fit_tatsuoka import load_literature
from .k4_crossmodel import REF, SIG_EV_PCT, SIG_SR, load_theta, metrics, run_dm04, run_norsand
from .sand import DATA_DIR, TIMS_P_A

DEPS_PCT = 0.05                                  # NorSand axial increment, % (K3's Weights.deps_pct)
POINT_END_PCT = 8.0                              # point tests without eps_peak run this far


def k3_curves():
    """The K3 curve set of fit_tatsuoka --lit, in its order (Tatsuoka 1986, then Lam & Tatsuoka, Fukushima)."""
    out = [DATA.load_curve(p) for p in sorted(glob.glob(os.path.join(DATA_DIR, "tatsuoka1986", "tats86_*.csv")))]
    return out + load_literature()


def vi_curves():
    return [DATA.load_curve(p) for p in sorted(glob.glob(os.path.join(DATA_DIR, "wang_toyoura", "wang34_Tokyo_*.csv")))]


def data_metrics(lc):
    """Peak metrics of a lab curve, by the k4_crossmodel rule on the sr branch; ev at the peak point interpolated from
    the eps_v branch (None outside it)."""
    ea, sr = np.abs(np.asarray(lc.eps_a_sr, float)), np.asarray(lc.sr, float)
    o = np.argsort(ea, kind="stable")
    ea, sr = ea[o], sr[o]
    if len(lc.eps_a_ev):
        xv, ev = np.abs(np.asarray(lc.eps_a_ev, float)), np.asarray(lc.eps_v, float)
        ov = np.argsort(xv, kind="stable")
        xv, ev = xv[ov], ev[ov]
        ev_on = np.where((ea >= xv[0]) & (ea <= xv[-1]), np.interp(ea, xv, ev), np.nan)
    else:
        ev_on = np.full_like(ea, np.nan)
    m = metrics("PS" if lc.kind.startswith("PS") else "TC", ea, sr, None, np.ones_like(sr), ev_on)
    m["ev_peak"] = None if not np.isfinite(m["ev_peak"]) else m["ev_peak"]
    m["ea_data_end"] = float(max(ea[-1], np.abs(lc.eps_a_ev).max() if len(lc.eps_a_ev) else 0.0))
    return m


def point_sr(phi_deg):
    s = math.sin(math.radians(phi_deg))
    return (1.0 + s) / (1.0 - s)


def specimens():
    """Every specimen: id, dataset, kind, sigma3, e, run length, and the data's peak metrics."""
    out = []
    for ds, curves in (("K3", k3_curves()), ("VI", vi_curves())):
        for lc in curves:
            d = data_metrics(lc)
            end = d["ea_data_end"]
            out.append(dict(id=lc.meta["id"], dataset=ds, kind="PS" if lc.kind.startswith("PS") else "TC",
                            p0=lc.sigma3, e=lc.e0, eps_a_end=end / 100.0, n=int(math.ceil(end / DEPS_PCT - 1e-9)),
                            weight=float(lc.meta.get("weight", 1.0)), data=d))
    pts = DATA.load_points(os.path.join(DATA_DIR, "tatsuoka1986", "point_tests.csv"), source="Tatsuoka 1986 Figs. 9 / 22")
    for k, pt in enumerate(pts):
        end = 2.0 * pt.eps_peak_pct + 1.0 if pt.eps_peak_pct is not None else POINT_END_PCT
        out.append(dict(id=f"tats86_pt{k:02d}_s{pt.sigma3:g}_e{pt.e:g}", dataset="K3", kind=pt.kind, p0=pt.sigma3, e=pt.e,
                        eps_a_end=end / 100.0, n=int(math.ceil(end / DEPS_PCT - 1e-9)), weight=1.0, point=True,
                        data=dict(sr_peak=point_sr(pt.phi_peak_deg), ea_peak=pt.eps_peak_pct, ev_peak=None)))
    return out


def judge(model, data):
    if not model:
        return dict(within=False, note="no complete model curve")
    d_sr = model["sr_peak"] - data["sr_peak"]
    d_ev = None if data.get("ev_peak") is None else model["ev_peak"] - data["ev_peak"]
    d_ea = None if data.get("ea_peak") is None else model["ea_peak"] - data["ea_peak"]
    ok = abs(d_sr) <= SIG_SR and (d_ev is None or abs(d_ev) <= SIG_EV_PCT)
    return dict(d_sr_peak=d_sr, d_ev_peak=d_ev, d_ea_peak=d_ea, within=bool(ok),
                sr_only=d_ev is None)


def _f(x, fmt="{:+.3f}"):
    return "—" if x is None else fmt.format(x)


def write_md(path, rep):
    L = [f"# K4, own data (O8: each model against its calibration data)", "",
         f"theta: `{rep['theta']}` ({rep['theta_source']}).", "",
         f"Gate: |d sr_peak| <= {SIG_SR:g} and |d ev_peak| <= {SIG_EV_PCT:g} % (K3's scatter model), on the diagonal "
         f"only: NorSand against the K3 set, DM04 `{REF}` against Verdugo & Ishihara (1996). p_a {TIMS_P_A:g} kPa.", ""]
    L += ["## Summary", "", "| model | data | specimens | within | sr only | median abs(d sr) | median abs(d ev) % | gated |",
          "|---|---|---|---|---|---|---|---|"]
    for row in rep["summary"]:
        L.append(f"| {row['model']} | {row['dataset']} | {row['n']} | {row['within']} | {row['sr_only']} | "
                 f"{_f(row['med_sr'], '{:.3f}')} | {_f(row['med_ev'], '{:.3f}')} | {'yes' if row['gated'] else 'no'} |")
    for ds, title in (("K3", "K3 set (NorSand's own data)"), ("VI", "Verdugo & Ishihara 1996 (DM04's own data)")):
        L += ["", f"## {title}", "",
              "| specimen | kind | sigma3 kPa | e | data sr / ea / ev | NorSand sr / ea / ev | d sr, d ev | DM04 sr / ea / ev | d sr, d ev |",
              "|---|---|---|---|---|---|---|---|---|"]
        for s in rep["specimens"]:
            if s["dataset"] != ds:
                continue
            cells = []
            for m in ("norsand", "dm04"):
                mm, j = s[m]["metrics"], s[m + "_judge"]
                if not mm:
                    cells += [f"{s[m]['status']}", "—"]
                    continue
                cells += [f"{mm['sr_peak']:.3f} / {mm['ea_peak']:.2f} / {mm['ev_peak']:+.3f}",
                          f"{_f(j['d_sr_peak'])}, {_f(j['d_ev_peak'])}{' ✓' if j['within'] else ''}"]
            d = s["data"]
            L.append(f"| {s['id']} | {s['kind']} | {s['p0']:g} | {s['e']:g} | {d['sr_peak']:.3f} / {_f(d['ea_peak'], '{:.2f}')} / "
                     f"{_f(d['ev_peak'])} | " + " | ".join(cells) + " |")
    bad = [f"{m} {s['id']}: {s[m]['status']}" for s in rep["specimens"] for m in ("norsand", "dm04")
           if s[m].get("status") not in ("ok",)]
    L += ["", "## Runs that did not finish cleanly", ""] + ([f"- {b}" for b in bad] or ["- none"])
    open(path, "w", encoding="utf-8").write("\n".join(L) + "\n")


def summarise(specs):
    rows = []
    for m, label in (("norsand", "NorSand"), ("dm04", f"DM04 {REF}")):
        for ds in ("K3", "VI"):
            js = [s[m + "_judge"] for s in specs if s["dataset"] == ds]
            dsr = [abs(j["d_sr_peak"]) for j in js if "d_sr_peak" in j]
            dev = [abs(j["d_ev_peak"]) for j in js if j.get("d_ev_peak") is not None]
            rows.append(dict(model=label, dataset=ds, n=len(js), within=sum(j["within"] for j in js),
                             sr_only=sum(bool(j.get("sr_only")) for j in js),
                             med_sr=float(np.median(dsr)) if dsr else None, med_ev=float(np.median(dev)) if dev else None,
                             gated=(m == "norsand") == (ds == "K3")))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--theta", required=True, help="a K3 result JSON or a JSON dict")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "out", "k4_own"))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", default=None, help="comma list of specimen ids (smoke)")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    theta, theta_src = load_theta(a.theta)
    specs = specimens()
    if a.only:
        keep = set(a.only.split(","))
        specs = [s for s in specs if s["id"] in keep]
    t0 = time.perf_counter()
    with ProcessPoolExecutor(a.workers) as ex:
        fn = ex.map(run_norsand, [(s, theta, s["n"]) for s in specs])
        fd = ex.map(run_dm04, [(s, REF, TIMS_P_A) for s in specs])
        ns = dict(fn)
        dm = {cid: r for (cid, _), r in fd}
    for s in specs:
        s["norsand"], s["dm04"] = ns[s["id"]], dm[s["id"]]
        s["norsand_judge"] = judge(s["norsand"]["metrics"], s["data"])
        s["dm04_judge"] = judge(s["dm04"]["metrics"], s["data"])
    rep = dict(label="K4 own data", theta=theta, theta_source=theta_src, ref=REF, p_a=TIMS_P_A,
               tolerance=dict(sr=SIG_SR, ev_pct=SIG_EV_PCT), specimens=specs, summary=summarise(specs),
               seconds=time.perf_counter() - t0)
    json.dump(rep, open(os.path.join(a.out, "k4_own.json"), "w", encoding="utf-8"), indent=1, default=float)
    write_md(os.path.join(a.out, "k4_own.md"), rep)
    print(f"K4 own data: {len(specs)} specimens in {rep['seconds']:.0f} s -> {a.out}")
    for row in rep["summary"]:
        print(f"  {row['model']:>16} on {row['dataset']}: {row['within']}/{row['n']} within"
              f"{' (gated)' if row['gated'] else ''}")


if __name__ == "__main__":
    main()
