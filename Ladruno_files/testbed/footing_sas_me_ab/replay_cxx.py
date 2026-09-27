"""WP-138: replay the footing's dumped Gauss-point states through the C++
`ladrunoSANISANDReplay` command (WP-127), one row = one committed state + the
REAL strain increment that point received next (footing_ab.py replay CSVs).

Runs under CPython 3.12 -S with FOOTING_BIN naming the dist/bin to test:

    FOOTING_BIN=<dir> python -S replay_cxx.py --csv <replay.csv> --out <json>
        [--scheme 1|129] [--extra "<material tokens>"] [--maxsub 2000]

Every row is replayed with -convention compressionPositive (the CSV's), -type
PlaneStrain, -dt dt_next, -prevIncrNorm prevIncrNorm, -primed 1. Output: one
JSON object per row (rc, census, returned state, rho_alpha of the returned
alpha, f before/after), plus -- when the CSV carries the committed next state
-- the replay-vs-analysis stress difference (a faithfulness check of the dump:
for the scheme that produced the run it must be ~0).
"""
import argparse
import csv
import json
import math
import os
import sys

_BIN = os.environ["FOOTING_BIN"]
os.add_dll_directory(_BIN)
sys.path.insert(0, _BIN)
sys.path.append(r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages")
import numpy as np            # noqa: E402
import opensees as ops        # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(os.path.abspath(_BIN)), ops.__file__

SAN = [264.32, 0.312885, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0,
       0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
MC, CC, NB, MM, E0, LAMC, XI, PATM = 1.3309, 0.71, 3.5, 0.005, 0.83, 0.027, 0.45, 101.0
STAT_NAMES = ["updates", "meCalls", "substeps", "accepted", "rejectedErr",
              "forcedAtDTmin", "forcedClampMc", "rejectedLowP", "abandonedLowP",
              "capHits", "entryPminClamps", "pnResets", "maxSubstepsOneUpdate",
              "lastSubsteps", "lastForcedAtDTmin", "lastAbandonedLowP", "lastCapHit"]


def nrm(v):
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2 + 2 * (v[3]**2 + v[4]**2 + v[5]**2))


def rho_alpha(alpha, e, p):
    an = nrm(alpha)
    if an == 0 or p <= 0:
        return float("nan")
    a = [x / an for x in alpha]
    M = np.array([[a[0], a[3], a[5]], [a[3], a[1], a[4]], [a[5], a[4], a[2]]])
    c3 = max(-1.0, min(1.0, math.sqrt(6.0) * float(np.trace(M @ M @ M))))
    g = 2 * CC / ((1 + CC) - (1 - CC) * c3)
    psi = e - (E0 - LAMC * (p / PATM) ** XI)
    return an / (math.sqrt(2.0 / 3.0) * (g * MC * math.exp(-NB * psi) - MM))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--scheme", type=int, default=1)
    ap.add_argument("--extra", default="")
    ap.add_argument("--maxsub", type=int, default=2000)
    ap.add_argument("--tolr", type=float, default=1.0e-7)
    a = ap.parse_args()
    extra = []
    for tok in a.extra.split():
        try:
            extra.append(int(tok))
        except ValueError:
            try:
                extra.append(float(tok))
            except ValueError:
                extra.append(tok)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("LadrunoSANISAND", 1, *SAN, a.scheme, 0, 1, 1.0e-7, a.tolr,
                   "-flipAlphaIn", "init", "-Pmin", 0.0101, "-maxSubsteps", a.maxsub,
                   "-Presidual", 0.0, "-honorTolR", 0, *extra)
    rows = list(csv.DictReader(open(a.csv, newline="")))
    out = []
    for r in rows:
        v = lambda k: [float(r[f"{k}_{i}"]) for i in range(6)]
        args = [1, "-convention", "compressionPositive",
                "-sigma", *v("sigma"), "-alpha", *v("alpha"), "-alphaIn", *v("alpha_in"),
                "-fabric", *v("z"), "-voidRatio", float(r["e"]),
                "-dStrain", *v("dStrain"), "-type", "PlaneStrain", "-trace", 0,
                "-dt", float(r["dt_next"]) if r.get("dt_next") else 1.0,
                "-primed", 1, "-prevIncrNorm",
                float(r["prevIncrNorm"]) if r.get("prevIncrNorm") else 0.0]
        res = ops.ladrunoSANISANDReplay(*args)
        res = list(res)
        ns = int(res[2])
        stats = dict(zip(STAT_NAMES, res[6:6 + ns]))
        b = 6 + ns
        sig = res[b:b + 6]; al = res[b + 6:b + 12]
        e_out, p_out, q_out, fb, fa, path, er = res[b + 24:b + 31]
        d = dict(element=int(r["element"]), gp=int(r["gp"]), step=int(r["step"]),
                 select=r["select"], rc=int(res[1]), stats=stats, sigma=sig,
                 alpha=al, e=e_out, p=p_out, q=q_out, f_before=fb, f_after=fa,
                 path=int(path), elastic_ratio=er,
                 rho_alpha_in=float(r["rho_alpha"]),
                 rho_alpha_out=rho_alpha(al, e_out, p_out))
        # WP-129 tail (IntScheme 129 binaries): 129, LSAS_COUNT, sasStats,
        # ratio before, ratio after, tangentEP(36). Trace is off (0 records).
        nrec, width = int(res[3]), int(res[4])
        t = b + 34 + nrec * width
        if len(res) > t + 1 and int(res[t]) == 129:
            nsas = int(res[t + 1])
            sas = res[t + 2:t + 2 + nsas]
            d["sas"] = sas
            d["sas_ratio_before"] = res[t + 2 + nsas]
            d["sas_ratio_after"] = res[t + 3 + nsas]
            d["sas_last_refuse_code"] = int(sas[27]) if nsas > 27 else None
            d["sas_substeps"] = int(sas[2])
        if r.get("sigma_next_0") and r["sigma_next_0"] != "nan":
            sn = v("sigma_next")
            ds = [s1 - s0 for s1, s0 in zip(sig, sn)]
            d["dsig_vs_analysis"] = nrm(ds)
            d["dsig_vs_analysis_rel"] = nrm(ds) / max(nrm(sn), 1e-12)
        out.append(d)
    with open(a.out, "w") as f:
        json.dump(dict(build=ops.ladrunoBuild().strip(), scheme=a.scheme,
                       extra=a.extra, csv=a.csv, rows=out), f, indent=1)
    rcs = [d["rc"] for d in out]
    print(f"{len(out)} rows, rc!=0: {sum(1 for x in rcs if x != 0)}, "
          f"max dsig vs analysis: "
          f"{max([d.get('dsig_vs_analysis', 0) for d in out], default=0):.3e}")


if __name__ == "__main__":
    main()
