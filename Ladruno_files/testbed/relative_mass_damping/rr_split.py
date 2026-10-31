"""Energy split (DW_rel, W_ref) and integrator guard checks for relativeMassDamping. Usage: py -3.12 rr_split.py <OpenSees.exe>"""
import subprocess
import sys
import pathlib

import numpy as np

EXE = sys.argv[1]
import tempfile
H = pathlib.Path(tempfile.gettempdir()) / "relative_mass_damping" / "rr_split_work"
H.mkdir(parents=True, exist_ok=True)
m, k, alpha, dt, n = 2.0, 800.0, 0.7, 0.005, 4000
t = np.arange(n + 1) * dt
vr = 0.3 * np.sin(2 * np.pi * 1.3 * t) + 0.1 * np.sin(2 * np.pi * 7.0 * t)
np.savetxt(H / "vr.txt", vr, fmt="%.17g")


def deck(name, integ="Newmark 0.5 0.25", ref=True, extra_load=False, split=True):
    d = H / name
    d.mkdir(exist_ok=True)
    (d / "vr.txt").write_bytes((H / "vr.txt").read_bytes())
    L = ["wipe", "model basic -ndm 1 -ndf 1", "node 1 0.0", "node 2 0.0 -mass 2.0", "fix 1 1",
         f"uniaxialMaterial Elastic 1 {k!r}", "element zeroLength 1 1 2 -mat 1 -dir 1",
         f"region 1 -nodeOnly 2 -rayleigh {alpha!r} 0.0 0.0 0.0",
         f"timeSeries Path 5 -dt {dt!r} -filePath vr.txt",
         # a base-like pulse so the structure is also loaded directly
         "timeSeries Path 6 -time {0.0 0.2 0.4 100.0} -values {0.0 50.0 0.0 0.0}", "pattern Plain 6 6 { load 2 1.0 }",
         f"recorder EnergyBalance -time -precision 15 -v2{' -relativeMassSplit' if split else ''} -file e2.txt",
         "recorder Node -file v.txt -time -precision 17 -node 2 -dof 1 vel"]
    if ref:
        L.append('puts "RR [relativeMassDamping 1 -vel 1 5 -node 2]"')
    if extra_load:
        L.append(f"pattern Plain 5 5 {{ load 2 {alpha * m!r} }}")
    L += ["constraints Plain", "numberer Plain", "system BandGeneral", "test NormDispIncr 1e-14 20",
          "algorithm Newton", f"integrator {integ}", "analysis Transient",
          f"set ok [analyze {n} {dt!r}]", 'puts "done $ok [getTime]"', "wipe"]
    (d / "main.tcl").write_text("\n".join(L) + "\n")
    r = subprocess.run([EXE, "main.tcl"], cwd=d, capture_output=True, text=True)
    return d, r.stdout + r.stderr, r.returncode


# --- split, recorder layout
d, txt, rc = deck("R_split")
hdr = [l for l in txt.splitlines() if "columns" in l]
print("R_split:", hdr, [l for l in txt.splitlines() if l.startswith(("RR", "done"))])
e = np.loadtxt(d / "e2.txt")
# time KE_ele KE_nod IE DW_ele DW_nod ULW DW_rel W_ref RES ERR
cols = ["time", "KE_ele", "KE_nod", "IE", "DW_ele", "DW_nod", "ULW", "DW_rel", "W_ref", "RES", "ERR"]
assert e.shape[1] == len(cols), e.shape
E = {c: e[:, i] for i, c in enumerate(cols)}
v = np.r_[0.0, np.loadtxt(d / "v.txt")[:, 1]]
rel = alpha * m * (v - vr) ** 2
ref = alpha * m * vr * (v - vr)
off_rel = np.r_[0.0, np.cumsum(0.5 * (rel[:-1] + rel[1:]) * dt)][1:]
off_ref = np.r_[0.0, np.cumsum(0.5 * (ref[:-1] + ref[1:]) * dt)][1:]
print(f"  DW_rel min {E['DW_rel'].min():.3e} (non-decreasing: {bool(np.all(np.diff(E['DW_rel']) >= -1e-15))}), end {E['DW_rel'][-1]:.6e}")
print(f"  W_ref end {E['W_ref'][-1]:.6e}; DW_nod end {E['DW_nod'][-1]:.6e}")
print(f"  DW_rel + W_ref - DW_nod: max {np.abs(E['DW_rel'] + E['W_ref'] - E['DW_nod']).max():.2e}")
print(f"  DW_rel vs offline trapezoid: {np.abs(E['DW_rel'] - off_rel).max() / np.abs(off_rel).max():.2e}; "
      f"W_ref vs offline: {np.abs(E['W_ref'] - off_ref).max() / np.abs(off_ref).max():.2e}")
print(f"  RES max {np.abs(E['RES']).max():.3e} (scale {np.abs(e[:, 1:9]).max():.3e})")

# --- byte identity of the -v2 layout without the flag and without a reference
d0, t0, _ = deck("A_nosplit", ref=False, extra_load=True, split=False)
d1, t1, _ = deck("A_split_flag", ref=False, extra_load=True, split=True)
a0 = np.loadtxt(d0 / "e2.txt"); a1 = np.loadtxt(d1 / "e2.txt")
print("A (load, no reference): columns without flag", a0.shape[1], "with flag", a1.shape[1],
      "; common columns equal:", bool(np.array_equal(a0[:, :7], a1[:, :7]) and np.array_equal(a0[:, -2:], a1[:, -2:])),
      "; split columns zero:", bool(np.all(a1[:, 7:9] == 0.0)))
E0 = a0
print(f"  A: DW_nod end {E0[-1, 5]:.6e}, ULW end {E0[-1, 6]:.6e}; R: DW_nod end {E['DW_nod'][-1]:.6e} ULW {E['ULW'][-1]:.6e}")
print(f"  A DW_nod - W_A(=ULW_A-ULW_R) = {E0[-1, 5] - (E0[-1, 6] - E['ULW'][-1]):.6e}  vs  R DW_nod {E['DW_nod'][-1]:.6e}")

# --- auto-enable: reference present at initialize, no flag
d2, t2, _ = deck("R_auto", split=False)
print("R_auto (no flag, reference defined before analyze):", [l for l in t2.splitlines() if "columns" in l])

# --- guard
for integ in ("HHT_TP 0.9", "Newmark 0.5 0.25", "HHT 0.9"):
    dg, tg, rc = deck("guard_" + integ.split()[0], integ=integ, split=False)
    lines = [l for l in tg.splitlines() if l.startswith(("FATAL", "done"))]
    print(f"guard {integ}: rc {rc}", lines[:2])
