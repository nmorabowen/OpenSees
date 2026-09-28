"""C++ SAS-ME (IntScheme 129) on the same refuser fan as the oracle, one replay per
trial; plus the raw replay vectors of a fixed set as a BYTE-IDENTITY baseline.

Run under the fork's CPython 3.12 with -S (no scipy there):
    py -3.12 -S cxx_fan.py <bin_dir> <out.json> [extra deck flags ...]
e.g. extra flags "-sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5" on a WP-151 build.

Prototype = the WP-138 E_B material: campaign set, IntScheme 129, TanType 0,
JacoType 1, TolF 1e-7, TolR 1e-4, -flipAlphaIn init, -Pmin 0.0101,
-maxSubsteps 2000, -Presidual 0, -honorTolR 0.
Jobs:
  fan    5 refuser states (data/refuser_states.csv) x 32 Fibonacci directions x
         {3e-6, 3e-5}, as a3_fan.py (Mandel-scaled plane-strain directions);
  ring   the 80 TIMs ring rows x {isoComp, shear} x {1e-5, 1e-4};
  repro  WP-128's smallest reproducer (sigma = 0.0101 I, dEps_yy in {1e-5, 1e-4, 3e-4}).
Output: every job's rc, sasStats tail (named), end state, and the full raw output
vector as float.hex (the byte-identity record)."""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.environ.get("R1_WORKTREE",
                   r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\sharp-chandrasekhar-d6ff72")
SAS_NAMES = [
    "updates", "elastic", "substeps", "accepted", "rejectedErr", "rejectedLowP",
    "rejectedNonPosH", "rejectedDrift", "rejectedAlpha", "elasticStages",
    "driftCorrections", "alphaInReseats", "hBrackets", "alphaProjected",
    "intersectFail", "refusals", "refStartF", "refStartAlpha", "refStartOther",
    "refDTmin", "refNonPosH", "refLowP", "refDrift", "refAlpha", "refCap",
    "maxSubstepsOneUpdate", "lastSubsteps", "lastRefuseCode", "maxAlphaRatio",
    "lastAlphaRatio", "lastF", "entryOverKappa", "rejectedReversal",
    "hFloored", "hSoftCapped", "reseatHeld",          # WP-151 (absent before)
]
CAMPAIGN = [264.32, 0.312885, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005,
            1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]


def fib_sphere(n):
    pts, ga = [], math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = 1.0 - 2.0 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1.0 - y * y))
        pts.append((r * math.cos(ga * i), y, r * math.sin(ga * i)))
    return pts


def deps_of(u, d):
    return [d * u[0], d * u[1], 0.0, d * u[2] * math.sqrt(2.0), 0.0, 0.0]


def load_refusers():
    rows = []
    for r in csv.DictReader(open(os.path.join(HERE, "data", "refuser_states.csv"), newline="")):
        g = lambda n: [float(r[f"{n}_{i}"]) for i in range(6)]
        rows.append(dict(leg=r["leg"], k=int(r["k"]), sigma=g("sigma"), alpha=g("alpha"),
                         alpha_in=g("alpha_in"), z=g("z"), e=float(r["e"])))
    return rows


def load_ring():
    rows = []
    for mesh in ("b8", "b16"):
        p = os.path.join(W, "Ladruno_implementation", "_tims_2d_model_requests_2026-09-25",
                         f"ring_points_{mesh}.csv")
        for r in csv.DictReader(open(p, newline="")):
            g = lambda n: [float(r[f"{n}_{i}"]) for i in range(6)]
            rows.append(dict(mesh=mesh, element=int(r["element"]), gp=int(r["gp"]),
                             sigma=g("sigma"), alpha=g("alpha"), alpha_in=g("alpha_in"),
                             z=g("z"), e=float(r["e"])))
    return rows


def main():
    bin_dir, out_path = sys.argv[1], sys.argv[2]
    extra = []
    for tok in sys.argv[3:]:
        try:
            extra.append(float(tok))
        except ValueError:
            extra.append(tok)
    site = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages"
    os.add_dll_directory(bin_dir)
    sys.path.insert(0, bin_dir)
    if site not in sys.path:
        sys.path.append(site)
    import opensees as ops
    here = os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__)))
    assert here == os.path.normcase(os.path.abspath(bin_dir)), ops.__file__
    ops.wipe()
    ops.nDMaterial("LadrunoSANISAND", 1, *CAMPAIGN, 129, 0, 1, 1.0e-7, 1.0e-4,
                   "-flipAlphaIn", "init", "-Pmin", 0.0101, "-maxSubsteps", 2000,
                   "-Presidual", 0.0, "-honorTolR", 0, *extra)

    def replay(st, de):
        args = [1, "-convention", "compressionPositive",
                "-sigma", *st["sigma"], "-alpha", *st["alpha"], "-alphaIn", *st["alpha_in"],
                "-fabric", *st["z"], "-voidRatio", st["e"], "-dStrain", *de,
                "-type", "3D", "-trace", 0, "-dt", 1.0, "-primed", 1, "-prevIncrNorm", 0.0]
        r = list(ops.ladrunoSANISANDReplay(*args))
        rc, nst, nrec, width = int(r[1]), int(r[2]), int(r[3]), int(r[4])
        i = 6 + nst
        stt = r[i:i + 34]
        i += 34 + nrec * width
        sas = None
        if len(r) > i + 1 and int(r[i]) == 129:
            n_sas = int(r[i + 1])
            sas = dict(zip(SAS_NAMES, r[i + 2:i + 2 + n_sas]))
        return dict(rc=rc, sas=sas, sigma=stt[0:6], alpha=stt[6:12], alpha_in=stt[12:18],
                    f_after=stt[28], hex=[float(x).hex() for x in r])

    out = dict(bin=bin_dir, extra=[str(x) for x in extra], fan=[], ring=[], repro=[])
    dirs = fib_sphere(32)
    for st in load_refusers():
        for d in (3e-6, 3e-5):
            for j, u in enumerate(dirs):
                o = replay(st, deps_of(u, d))
                o.update(leg=st["leg"], k=st["k"], delta=d, idir=j)
                out["fan"].append(o)
    for row in load_ring():
        for d in (1e-5, 1e-4):
            for pn, de in (("isoComp", [d, d, 0, 0, 0, 0]), ("shear", [0, 0, 0, d, 0, 0])):
                o = replay(row, de)
                o.update(mesh=row["mesh"], element=row["element"], gp=row["gp"], probe=pn, delta=d)
                out["ring"].append(o)
    z6 = [0.0] * 6
    for d in (1e-5, 1e-4, 3e-4):
        st = dict(sigma=[0.0101, 0.0101, 0.0101, 0, 0, 0], alpha=z6, alpha_in=z6, z=z6, e=0.6944)
        o = replay(st, [0, d, 0, 0, 0, 0])
        o.update(delta=d)
        out["repro"].append(o)
    json.dump(out, open(out_path, "w"))
    nref = sum(1 for o in out["fan"] if o["rc"] != 0)
    print("fan refusals", nref, "/", len(out["fan"]), "| ring rc!=0",
          sum(1 for o in out["ring"] if o["rc"] != 0), "/", len(out["ring"]))


if __name__ == "__main__":
    main()
