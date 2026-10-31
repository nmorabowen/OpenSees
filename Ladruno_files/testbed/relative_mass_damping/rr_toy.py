"""Check of `relativeMassDamping` on the drag toy.

Runs per (pair, record), all from drag_toy.deck():
  ref, abs, A_ff, A_rec          old binary (dist/bin, 1652f945c) and new binary (worktree, pin + feature)
  R_ff, R_rec                    new binary: A deck with the compensating pattern replaced by
                                 `relativeMassDamping 1 -vel 1 13 -region 1` (same series 13)
Reports: R vs ref, R vs A (oracle), A(new) vs A(old) bitwise, ref/abs(new) vs (old) bitwise, energy closure.
"""
import hashlib
import pathlib
import re
import subprocess
import sys

import numpy as np

import os
import tempfile
# drag_toy.py lives in the Paper 1 repository (engstruct-2026/20_models/ssi/drag)
DRAG = pathlib.Path(os.environ.get("DRAG_TOY_DIR", r"C:\Users\nmb\Documents\gitAPE\Epistemic Uncertanty\papers\engineering-structures\review\engstruct-2026\20_models\ssi\drag"))
sys.path.insert(0, str(DRAG))
import drag_toy as T  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[3]
OLD = pathlib.Path(os.environ.get("RMD_OLD_EXE", r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin\OpenSees.exe"))   # a binary without the feature
NEW = pathlib.Path(os.environ.get("RMD_NEW_EXE", str(REPO / "dist" / "bin" / "OpenSees.exe")))
WORK = pathlib.Path(os.environ.get("RMD_WORK", str(pathlib.Path(tempfile.gettempdir()) / "relative_mass_damping" / "rr_work")))


def run(exe, d):
    r = subprocess.run([str(exe), "main.tcl"], cwd=d, capture_output=True, text=True)
    txt = r.stdout + r.stderr
    if "done 0" not in txt:
        raise RuntimeError(f"{d}: {txt[-2000:]}")
    return txt


def md5(p):
    return hashlib.md5(p.read_bytes()).hexdigest()


def make(run_name, pair, k, d, energy=True):
    base = {"R_ff": "A_ff", "R_rec": "A_rec"}.get(run_name, run_name)
    out = T.deck(base, pair, k, d)
    s = (d / "main.tcl").read_text(encoding="utf-8")
    if run_name.startswith("R_"):
        s, n = re.subn(r"pattern Plain 13 13 \{\n(?:  load .*\n)+\}\n", "set rr [relativeMassDamping 1 -vel 1 13 -region 1]\nputs \"RR $rr\"\n", s)
        assert n == 1
    if energy:
        s = s.replace("constraints Penalty", "recorder EnergyBalance -time -file energy.txt\nconstraints Penalty", 1)
    (d / "main.tcl").write_text(s, encoding="utf-8")
    return out


def energy_res(d):
    e = np.loadtxt(d / "energy.txt")
    # columns: time KE IE DW ULW RES ERR
    ke, ie, dw, ulw, res = e[:, 1], e[:, 2], e[:, 3], e[:, 4], e[:, 5]
    scale = np.abs(np.c_[ke, ie, dw, ulw]).max()
    return float(np.abs(res).max() / scale), float(dw[-1]), float(ulw[-1])


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    print("old", OLD, "\nnew", NEW)
    for sta in T.STATIONS:
        t, acc = T.C.read_record(sta, 1)
        ag = np.zeros(T.NSTEP + 1); ag[:len(acc)] = acc[:T.NSTEP + 1, 0]
        vrec = T.C.cumtrapz0(ag[:, None])[:, 0]
        urec = T.C.cumtrapz0(vrec[:, None])[:, 0]
        dff = WORK / sta / "ff"
        dff.mkdir(parents=True, exist_ok=True)
        T.write_series(dff / "vrec.txt", vrec)
        T.ff_deck(dff)
        run(OLD, dff)
        ff = {q: np.loadtxt(dff / f"{q}.txt") for q in ("disp", "vel", "accel")}
        tt = ff["disp"][:, 0]
        uff, vff, aff = ff["disp"][:, 1], ff["vel"][:, 1], ff["accel"][:, 1]
        for case, (pair, T1) in T.CASES.items():
            k = T.storey_k(T1)
            res, dirs = {}, {}
            jobs = [("ref", OLD), ("abs", OLD), ("A_ff", OLD), ("A_rec", OLD),
                    ("ref", NEW), ("abs", NEW), ("A_ff", NEW), ("A_rec", NEW), ("R_ff", NEW), ("R_rec", NEW)]
            for rn, exe in jobs:
                tagx = f"{rn}_{'new' if exe == NEW else 'old'}"
                d = WORK / sta / case / tagx
                d.mkdir(parents=True, exist_ok=True)
                T.write_series(d / "vrec.txt", vrec)
                T.write_series(d / "urec.txt", urec)
                T.write_series(d / "vff.txt", T.series_on_grid(tt, vff))
                T.write_series(d / "uff.txt", T.series_on_grid(tt, uff))
                T.write_series(d / "aff.txt", T.series_on_grid(tt, aff))
                out = make(rn, pair, k, d)
                txt = run(exe, d)
                if rn.startswith("R_"):
                    rr = [l for l in txt.splitlines() if l.startswith("RR ")]
                    print(f"  {sta} {case} {tagx}: {rr}  (alpha*sum m = {pair[0] * (T.M_FOUND + 3 * T.M_FLOOR):.6g})")
                res[tagx] = T.quantities(T.load(d, out), rn if not rn.startswith("R_") else "A_ff", pair, k, uff, vff)
                dirs[tagx] = d
            r = res["ref_old"]
            for x in ("ref_new", "abs_new", "A_ff_new", "A_rec_new"):
                same = md5(dirs[x] / "disp.txt") == md5(dirs[x.replace("_new", "_old")] / "disp.txt")
                same_v = md5(dirs[x] / "vel.txt") == md5(dirs[x.replace("_new", "_old")] / "vel.txt")
                same_e = md5(dirs[x] / "energy.txt") == md5(dirs[x.replace("_new", "_old")] / "energy.txt")
                print(f"  {sta} {case} {x:10s} bitwise vs old binary: disp {same} vel {same_v} energy {same_e}")
            for x in ("A_ff_old", "A_rec_old", "R_ff_new", "R_rec_new", "abs_old"):
                q = res[x]
                print(f"  {sta} {case} {x:10s} vs ref: roof {T.err(q['roof'], r['roof']):.2e} V {T.err(q['V'], r['V']):.2e} "
                      f"found {T.err(q['found'], r['found']):.2e}")
            for xr, xa in (("R_ff_new", "A_ff_old"), ("R_rec_new", "A_rec_old")):
                qa, qr = res[xa], res[xr]
                print(f"  {sta} {case} {xr:10s} vs {xa}: roof {T.err(qr['roof'], qa['roof']):.2e} V {T.err(qr['V'], qa['V']):.2e} "
                      f"found {T.err(qr['found'], qa['found']):.2e}")
            for x in ("abs_old", "A_ff_old", "A_rec_old", "R_ff_new", "R_rec_new"):
                rel, dw, ulw = energy_res(dirs[x])
                print(f"  {sta} {case} {x:10s} energy: max|RES|/max component {rel:.2e}  DW_end {dw:.6g}  ULW_end {ulw:.6g}")


if __name__ == "__main__" and "--sdof" not in sys.argv:
    main()


def sdof():
    """T1: SDOF m u'' + alpha m (u' - v_ref) + k u = 0, Newmark average acceleration, against the closed recursion."""
    m, k, alpha, dt, n = 2.0, 800.0, 0.7, 0.005, 4000
    t = np.arange(n + 1) * dt
    vr = 0.3 * np.sin(2 * np.pi * 1.3 * t) + 0.1 * np.sin(2 * np.pi * 7.0 * t)
    d = WORK / "sdof"
    d.mkdir(parents=True, exist_ok=True)
    np.savetxt(d / "vr.txt", vr, fmt="%.17g")
    for mode in ("R", "A"):
        L = ["wipe", "model basic -ndm 1 -ndf 1", "node 1 0.0", "node 2 0.0 -mass 2.0", "fix 1 1",
             f"uniaxialMaterial Elastic 1 {k!r}", "element zeroLength 1 1 2 -mat 1 -dir 1",
             f"region 1 -nodeOnly 2 -rayleigh {alpha!r} 0.0 0.0 0.0",
             f"timeSeries Path 5 -dt {dt!r} -filePath vr.txt"]
        L += ["puts \"RR [relativeMassDamping 1 -vel 1 5 -node 2]\""] if mode == "R" else [f"pattern Plain 5 5 {{ load 2 {alpha * m!r} }}"]
        L += [f"recorder Node -file {mode}_disp.txt -time -precision 17 -node 2 -dof 1 disp",
              f"recorder Node -file {mode}_rf.txt -time -precision 17 -node 2 -dof 1 rayleighForces",
              f"recorder Node -file {mode}_nrf.txt -time -precision 17 -node 2 -dof 1 nodalRayleighForces",
              "recorder EnergyBalance -time -file " + f"{mode}_energy.txt",
              "constraints Plain", "numberer Plain", "system BandGeneral", "test NormDispIncr 1e-14 10",
              "algorithm Linear", "integrator Newmark 0.5 0.25", "analysis Transient",
              f"set ok [analyze {n} {dt!r}]", 'puts "done $ok [getTime]"', "wipe"]
        (d / "main.tcl").write_text("\n".join(L) + "\n", encoding="utf-8")
        print("  sdof", mode, [l for l in run(NEW, d).splitlines() if l.startswith(("RR", "done"))])
    # closed recursion
    u = np.zeros(n + 1); v = np.zeros(n + 1); a = np.zeros(n + 1)
    b, g = 0.25, 0.5
    keff = k + alpha * m * g / (b * dt) + m / (b * dt * dt)
    for i in range(n):
        rhs = (alpha * m * vr[i + 1]
               + m * (u[i] / (b * dt * dt) + v[i] / (b * dt) + (1 / (2 * b) - 1) * a[i])
               + alpha * m * (g / (b * dt) * u[i] + (g / b - 1) * v[i] + dt * (g / (2 * b) - 1) * a[i]))
        u[i + 1] = rhs / keff
        a[i + 1] = (u[i + 1] - u[i]) / (b * dt * dt) - v[i] / (b * dt) - (1 / (2 * b) - 1) * a[i]
        v[i + 1] = v[i] + dt * ((1 - g) * a[i] + g * a[i + 1])
    R = np.loadtxt(d / "R_disp.txt")[:, 1]; A = np.loadtxt(d / "A_disp.txt")[:, 1]
    rf = np.loadtxt(d / "R_rf.txt")[:, 1]
    print(f"  sdof R vs recursion {np.abs(R - u[1:]).max() / np.abs(u).max():.2e}; R vs A {np.abs(R - A).max() / np.abs(A).max():.2e}")
    nrf = np.loadtxt(d / "R_nrf.txt")[:, 1]
    rf_exp = alpha * m * (v[1:] - vr[1:])
    print(f"  sdof rayleighForces (reaction flag 2) vs +alpha m (v - v_ref) {np.abs(rf - rf_exp).max() / np.abs(rf_exp).max():.2e}")
    print(f"  sdof nodalRayleighForces vs -alpha m (v - v_ref) {np.abs(nrf + rf_exp).max() / np.abs(rf_exp).max():.2e}")
    for mode in ("R", "A"):
        e = np.loadtxt(d / f"{mode}_energy.txt")
        print(f"  sdof {mode} energy max|RES|/max comp {np.abs(e[:, 5]).max() / np.abs(e[:, 1:5]).max():.2e}  DW_end {e[-1, 3]:.6g} ULW_end {e[-1, 4]:.6g}")


if __name__ == "__main__" and "--sdof" in sys.argv:
    sdof()
