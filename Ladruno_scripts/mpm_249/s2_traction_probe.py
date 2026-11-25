"""#1110 probe: S2's column with the surcharge carried by `mpm traction` (232l, P10b).

python3.12 Ladruno_scripts/mpm_249/s2_traction_probe.py <bindir> <repo> <mat> <basis> [nsur]

The grid's top face is put AT the soil top (y = l0), because `mpm traction` loads the nodes of a
grid boundary face row. Gravity is applied as the unscaled per-particle body force (B11), so the
load factor lambda drives only the traction: 0 through stage 0, the flip and the hold, then a
ramp to 1 over nsur steps, as S2's body-force ramp. Everything else is s2_run's timeline.
"""
import json
import os
import sys

bindir, repo, mat, basis = sys.argv[1:5]
nsur = int(sys.argv[5]) if len(sys.argv) > 5 else 10
os.add_dll_directory(bindir)
sys.path.insert(0, bindir)
sys.path.insert(0, os.path.join(repo, "tests"))
import opensees as ops  # noqa: E402
import ladruno_mpm_soil_s_harness as H  # noqa: E402

print("BUILD", ops.ladrunoBuild(), flush=True)
rec = H.load("s2_column.json")
leg = rec["legs"][f"{mat}_nsur{nsur}"]
l0, hcell, G, RHO = 10.0, 0.5, H.G, H.RHO
q = float(leg.get("surcharge", RHO * G * l0))
rgl = RHO * G * l0
DT = 0.1
t_on = DT * 12          # 10 gravity steps + 2 hold steps
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.nDMaterial(*H.deck_args(rec["decks"][mat]["cmd"]))
ops.mpm("grid", H.HOST, "-origin", 0.0, 0.0, "-h", hcell, "-cells", 2, int(round(l0 / hcell)),
        "-basis", basis)
ops.mpm("mode", "quasiStatic")
ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, l0, "-ppc", 2, "-mat", 1, "-rho", RHO)
ops.mpm("fix", "xmin", "x")
ops.mpm("fix", "xmax", "x")
ops.mpm("fix", "ymin", "y")
ops.mpm("gravity", 0.0, 0.0)
ops.mpm("traction", "sur", "ymax", 0.0, l0, 1.0, l0, 0.0, -q)
ops.timeSeries("Path", H.TS_TAG, "-time", 0.0, t_on, t_on + nsur * DT, 1.0e9,
               "-values", 0.0, 0.0, 1.0, 1.0)
ops.mpm("timeSeries", H.TS_TAG)
ops.constraints("Plain")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormDispIncr", 1.0e-10, 25)
ops.algorithm("Newton")
ops.integrator("LoadControl", DT)
ops.analysis("Static")
ops.mpm("stage", 1, 0)

x0 = H._l(ops.mpm("particleData", "x"))
N = len(x0) // 2
Y = [x0[2 * k + 1] for k in range(N)]
out = {"mat": mat, "basis": basis, "nsur": nsur, "carrier": "mpm traction ymax (grid top = soil top)",
       "build": ops.ladrunoBuild(), "checkpoints": {}}


def check(name):
    st = H._column_state(ops)
    ey, ex = H._column_errors(st, Y, leg[name]["particles"], rgl)
    out["checkpoints"][name] = {"syy_max": ey, "sxx_rms": ex}
    print(name, out["checkpoints"][name], flush=True)


def step(lam_g):
    ops.mpm("particleBodyForce", *([0.0, -G * lam_g] * N))
    ops.mpm("analyze", 1, DT)


for i in range(1, 11):
    step(i / 10.0)
check("end_stage0")
ops.mpm("stage", 1, 1)
check("after_flip")
step(1.0)
step(1.0)
check("after_hold")
for i in range(1, nsur + 1):
    step(1.0)
check("after_surcharge")
# the top particle row and the row-wise worst error after the surcharge
st = H._column_state(ops)
rows = leg["after_surcharge"]["particles"]
ytop = max(Y)
worst = {}
for k in range(N):
    r = min(rows, key=lambda qq: abs(qq["Y"] - Y[k]))
    e = (st["sig"][3 * k + 1] - r["syy"]) / rgl
    yk = round(Y[k], 4)
    worst.setdefault(yk, []).append(e)
out["after_surcharge_syy_signed_by_row"] = {str(k): sum(worst[k]) / len(worst[k]) for k in sorted(worst, reverse=True)}
out["passed"] = all(c["syy_max"] <= 0.03 and c["sxx_rms"] <= 0.02
                    for n, c in out["checkpoints"].items() if n != "after_hold")
print("RESULT", json.dumps(out), flush=True)
