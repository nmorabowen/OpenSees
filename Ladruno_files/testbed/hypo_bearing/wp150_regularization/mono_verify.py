"""WP-150: verify a fitted monotonic set (out_mono_<label>.json) on the C++ LadrunoSANISAND (SAS-ME) with R1 OFF and ON,
and its cyclic undrained sanity (CTXu at e 0.6944, p0 100, CSR 0.15 / 0.20).
    python mono_verify.py <label> <dist/bin>
Writes out_mono_<label>_verify.md.  The C++ runs in a CPython 3.12 -S subprocess (gate0_cxx_driver.py)."""
import json, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PY = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe"
SITE = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages"
label, bin_dir = sys.argv[1], sys.argv[2]
fit = json.load(open(os.path.join(HERE, f"out_mono_{label}.json")))
nb, A0, nd, h0 = fit["pars"]
params = [264.32, 0.312885, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005, h0, 0.968, nb, A0, nd, 12.5, 1100.0, 2.0]
R1 = ["-sasHFloor", 1.0, "-sasReseatHyst", 1.0, "-sasSoftCap", 0.5]
tests = [[f"{k}_p{int(p)}", {"TX": "TXd", "PS": "PSd"}[k], 0.6944, float(p), 0.25] for k in ("TX", "PS") for p in (10, 50, 150, 500)]
tests += [["TXu_comp_p100", "TXu", 0.6944, 100.0, 0.10], ["TXu_ext_p100", "TXu", 0.6944, 100.0, -0.10]]
tests += [["CTXu_csr0.20", "CTXu", 0.6944, 100.0, 0.20], ["CTXu_csr0.15", "CTXu", 0.6944, 100.0, 0.15]]
spec = dict(params=params, tolr=1e-4, nstep=1000, ncyc=40, variants={"cxx_off": [], "cxx_on": R1}, tests=tests)
sp = os.path.join(HERE, f"mono_{label}_cxx_spec.json"); op = os.path.join(HERE, f"out_mono_{label}_cxx.json")
json.dump(spec, open(sp, "w"), indent=1)
env = dict(os.environ); env.pop("PYTHONPATH", None); env["PYTHONIOENCODING"] = "utf-8"
cp = subprocess.run([PY, "-S", "-u", os.path.join(HERE, "gate0_cxx_driver.py"), bin_dir, SITE, sp, op],
                    capture_output=True, text=True, env=env)
if cp.returncode != 0:
    raise SystemExit(cp.stdout[-3000:] + cp.stderr[-3000:])
c = json.load(open(op))
# FULL-resolution exact oracle, compared AT the C++ output points (a coarse-grid linear interpolant of the C++ in the
# steep first 0.1 % of strain otherwise reads as a spurious % error; the curves agree to <= 0.2 % at matched strain).
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", "Ladruno_scripts")))
from dataclasses import replace
from sanisand_reference import CAMPAIGN
from sanisand_reference.driver import isotropic_state
from sanisand_reference.integrator import Control, integrate, path_table
from sanisand_reference.ring import ring_variants
Pm = replace(CAMPAIGN, e_init=0.6944, nb=nb, A0=A0, nd=nd, h0=h0)
Om = ring_variants()["uw_model"]
orc = {}
for kind in ("TX", "PS"):
    for p0 in (10.0, 50.0, 150.0, 500.0):
        mask = (True, False, False, True, True, True) if kind == "TX" else (True, False, True, True, True, True)
        res = integrate(isotropic_state(p0, 0.6944), Control(mask, (0.25, 0, 0, 0, 0, 0)), Pm, Om, rtol=1e-9)
        tab = path_table(res, Pm, Om)
        orc[f"{kind}_p{int(p0)}"] = dict(eps_a=[t["eps"][0] for t in tab], q=[t["q"] for t in tab])


def dq(x, y):
    # x = reference (dense), y = compared (evaluated at its own points)
    ey, qy = np.array(y["eps_a"]), np.array(y["q"])
    m = ey <= min(x["eps_a"][-1], ey[-1])
    return float(np.max(np.abs(np.interp(ey[m], x["eps_a"], x["q"]) - qy[m])) / max(x["q"]))


L = [f"## C++ verification of `{label}` (nb {nb:.3f}, A0 {A0:.3f}, nd {nd:.3f}, h0 {h0:.3f}); build {c['build']}", "",
     "| test | max \\|Δq\\|/q_max, C++ R1 off vs exact oracle | R1 on vs off |", "|---|---|---|"]
w1 = w2 = 0.0
for t in orc:
    a, b = c["cxx_off"][t], c["cxx_on"][t]
    d1, d2 = dq(orc[t], a), dq(a, b); w1, w2 = max(w1, d1), max(w2, d2)
    L.append(f"| {t} | {d1:.2e} | {d2:.2e} |")
L += ["", f"worst: C++ vs oracle {w1:.2e}; R1 on vs off {w2:.2e}", "",
      "| TXu monotonic, e 0.6944, p0 100 (R1 on) | p′_min (kPa) at ε_a | p′ at \|ε_a\| 1 / 5 % |", "|---|---|---|"]
for t in ("TXu_comp_p100", "TXu_ext_p100"):
    r = c["cxx_on"][t]; pp = np.array(r["p"]); ee = np.abs(np.array(r["eps_a"])); i = int(np.argmin(pp))
    L.append(f"| {t} | {pp[i]:.1f} at {ee[i]:.4f} | {np.interp(0.01, ee, pp):.0f} / {np.interp(0.05, ee, pp):.0f} |")
L += ["",
      "| CTXu (e 0.6944, p0 100) | N at 5 % DA, R1 off | R1 on | p′_min off / on |", "|---|---|---|---|"]
for t in ("CTXu_csr0.20", "CTXu_csr0.15"):
    a, b = c["cxx_off"][t], c["cxx_on"][t]
    L.append(f"| {t} | {a['N']} | {b['N']} | {a['p_min']:.1f} / {b['p_min']:.1f} |")
txt = "\n".join(L) + "\n"
print(txt)
open(os.path.join(HERE, f"out_mono_{label}_verify.md"), "w", encoding="utf-8").write(txt)
