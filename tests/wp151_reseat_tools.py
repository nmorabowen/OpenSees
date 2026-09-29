"""WP-151 helpers (not collected: no `test_` prefix): the replay set for the R1
tests and for its byte-identity baseline, so both use EXACTLY the same jobs and
row format.

Material = the WP-138 E_B deck's LadrunoSANISAND: the TIMs campaign set,
IntScheme 129 (SAS-ME), TanType 0, JacoType 1, TolF 1e-7, TolR 1e-4,
-flipAlphaIn init, -Pmin 0.0101, -maxSubsteps 2000, -Presidual 0, -honorTolR 0.

Jobs (all COMPRESSION POSITIVE, tensor components, engineering-shear strain):
  fan    the five loadingNonPosH refusers at the Esmeralda walls
         (tests/data/wp151_refuser_states.csv) x 32 Fibonacci directions in the
         plane-strain strain space (Mandel-scaled) x {3e-6, 3e-5};
  ring   the 80 TIMs ring rows x {isoComp, shear} x {1e-5, 1e-4};
  repro  WP-128's smallest reproducer, dEps_yy in {1e-5, 1e-4, 3e-4}.
A ROW is [rc, then float.hex of sigma(6) alpha(6) alpha_in(6) z(6) e, and the
sasStats substeps, lastRefuseCode, alphaInReseats, rejectedReversal].
"""
import csv
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.join(os.path.dirname(_HERE), "Ladruno_scripts")
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

import sanisand_replay as sr  # noqa: E402

DATA = os.path.join(_HERE, "data")
REFUSERS_CSV = os.path.join(DATA, "wp151_refuser_states.csv")
BASELINE = os.path.join(DATA, "wp151_sasme_byteid_baseline.json")
ORACLE_FAN = os.path.join(DATA, "wp151_oracle_fan.json")

P = list(sr.CAMPAIGN_PARAMS)
EB_OPTS = (129, 0, 1, 1.0e-7, 1.0e-4, "-flipAlphaIn", "init", "-Pmin", 0.0101,
           "-maxSubsteps", 2000, "-Presidual", 0.0, "-honorTolR", 0)
CONE = math.sqrt(2.0 / 3.0) * P[9]          # yield-cone radius sqrt(2/3) m

# prototype tag -> extra deck flags
R1_ON = ("-sasHFloor", 1.0, "-sasReseatHyst", 1.0, "-sasSoftCap", 0.5)
PROTOS = {
    11: (),                                                        # DM04 (flags absent)
    12: ("-sasHFloor", 0.0, "-sasReseatHyst", 0.0, "-sasSoftCap", 0.0),   # flags given as 0
    13: R1_ON,                                                     # R1 as recommended
    14: ("-sasHFloor", 1.0),                                       # floor alone
    15: ("-sasReseatHyst", 1.0),                                   # hysteresis alone
    16: ("-sasHFloor", 1.0, "-sasReseatHyst", 1.0),                # floor + hysteresis, no cap
}


def define(ops, tags=None, protos=PROTOS):
    for tag in (tags or protos):
        ops.nDMaterial("LadrunoSANISAND", tag, *P, *EB_OPTS, *protos[tag])


def fib_sphere(n):
    pts, ga = [], math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = 1.0 - 2.0 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1.0 - y * y))
        pts.append((r * math.cos(ga * i), y, r * math.sin(ga * i)))
    return pts


def deps_of(u, d):
    return [d * u[0], d * u[1], 0.0, d * u[2] * math.sqrt(2.0), 0.0, 0.0]


def refusers():
    rows = []
    with open(REFUSERS_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            g = lambda n, r=r: [float(r[f"{n}_{i}"]) for i in range(6)]
            rows.append(dict(leg=r["leg"], k=int(r["k"]), element=int(r["element"]),
                             gp=int(r["gp"]), sigma=g("sigma"), alpha=g("alpha"),
                             alpha_in=g("alpha_in"), z=g("z"), e=float(r["e"]),
                             deps_last=g("deps")))
    return rows


def jobs():
    """[(kind, key, state, dstrain)] in a fixed order."""
    out = []
    dirs = fib_sphere(32)
    for st in refusers():
        for d in (3e-6, 3e-5):
            for j, u in enumerate(dirs):
                out.append(("fan", f"{st['leg']}/{st['k']}/{d:g}/{j}", st, deps_of(u, d)))
    for path in sr.RING_CSVS:
        mesh = "b8" if "b8" in os.path.basename(path) else "b16"
        for row in sr.read_ring_csv(path):
            for d in (1e-5, 1e-4):
                for pn, de in sr.probes(d).items():
                    out.append(("ring", f"{mesh}/{row['element']}/{row['gp']}/{pn}/{d:g}", row, de))
    z6 = [0.0] * 6
    for d in (1e-5, 1e-4, 3e-4):
        st = dict(sigma=[0.0101, 0.0101, 0.0101, 0.0, 0.0, 0.0], alpha=z6, alpha_in=z6, z=z6,
                  e=P[2])
        out.append(("repro", f"repro/{d:g}", st, [0.0, d, 0.0, 0.0, 0.0, 0.0]))
    return out


def replay(ops, tag, st, de):
    return sr.replay(ops, tag, st["sigma"], st["alpha"], st["alpha_in"], st["z"], st["e"], de,
                     "compressionPositive", trace=0)


def row(o):
    sas = o["sas"] or {}
    vals = (list(o["sigma"]) + list(o["alpha"]) + list(o["alpha_in"]) + list(o["z"]) + [o["e"]]
            + [sas.get("substeps", float("nan")), sas.get("lastRefuseCode", float("nan")),
               sas.get("alphaInReseats", float("nan")), sas.get("rejectedReversal", float("nan"))])
    return [int(o["rc"])] + [float(v).hex() for v in vals]


def run_rows(ops, tag):
    return {key: row(replay(ops, tag, st, de)) for (_k, key, st, de) in jobs()}


def rows_equal(cur, ref):
    """EXACT on win32 (the baseline's platform, MSVC); elsewhere the fork's 1e-6
    cross-platform floor on floats (as test_ladruno_sanisand_sasme.py), rc exact."""
    if sys.platform == "win32":
        return cur == ref
    if len(cur) != len(ref) or cur[0] != ref[0]:
        return False
    xs = [float.fromhex(x) for x in cur[1:]]
    ys = [float.fromhex(y) for y in ref[1:]]
    for x, y in zip(xs, ys):
        if math.isnan(x) and math.isnan(y):
            continue
        if abs(x - y) > 1e-6 * max(abs(y), 1.0):
            return False
    return True
