"""K4 (plan 144 §5.2), PM4Sand leg: single-element plane-strain compression in OpenSees, drained and undrained.

Runs under the interpreter the build targets (CPython 3.12, BUILD_GOTCHAS §4), which has no SciPy; the NorSand and
DM04 legs run under 3.11 in k4_crossmodel.py, which reads the JSON written here.

PM4Sand (Boulanger & Ziotopoulou) is a 2D plane-strain model with its OWN critical state line (Bolton's relative
dilatancy index, Q and R): it cannot share the fork CSL. It is matched to DM04's Toyoura set where it can be:
  Dr     = (e_max - e) / (e_max - e_min), Toyoura e_max 0.977, e_min 0.597 (Ishihara 1993);
  G0     = DM04 G0 (2.97 - e)^2 / (1 + e)  (PM4Sand: G = G0 p_a sqrt(p/p_a));
  phi_cv = asin(3 M / (6 + M)) with DM04 M = 1.25 (31.1 deg);  nu = DM04 0.05;  p_a = 101 kPa (TIMs).
  hpo    = 0.4 (contraction rate; it matters mostly for cyclic loading, monotonic paths are insensitive to it).
So PM4Sand is a secondary, mechanism-level comparison; DM04 shares the CSL exactly.

One SSPquad, unit square, thickness 1, plane strain. Axial = y (compression), lateral = x. Stage 0 (elastic)
consolidation to sigma_xx = sigma_yy = -p0, then updateMaterialStage 1.
  PS  : lateral force held (loadConst), top displaced (DisplacementControl), equalDOF keeps faces plane.
  PSU : constant volume, eps_xx = -eps_yy, both faces displaced by Path series from the consolidated state.
Output stress is PM4Sand's in-plane (sigma_xx, sigma_yy, tau_xy); out-of-plane sigma_zz is not a PM4Sand output,
so the comparison uses in-plane quantities: s = (s1 + s3)/2, t = (s1 - s3)/2, s1/s3, eps_v (= in-plane, eps_zz = 0).

    py -3.12 -S k4_pm4sand.py --dist <worktree>/dist/bin --cases cases.json --out pm4sand.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

E_MAX, E_MIN = 0.977, 0.597


def pm4_params(e, sand):
    Dr = (E_MAX - e) / (E_MAX - E_MIN)
    G0 = sand["G0"] * (2.97 - e) ** 2 / (1.0 + e)
    phi_cv = math.degrees(math.asin(3.0 * sand["M"] / (6.0 + sand["M"])))
    return dict(Dr=Dr, G0=G0, hpo=0.4, phi_cv=phi_cv, nu=sand["nu"], p_a=sand["p_a"])


def run_case(ops, case, sand, n=400):
    kind, p0, e, eps_end = case["kind"], case["p0"], case["e"], case["eps_a_end"]
    pr = pm4_params(e, sand)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], start=1):
        ops.node(tag, float(x), float(y))
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(4, 1, 0)
    ops.equalDOF(3, 2, 1)       # right face stays vertical (node 3 retained)
    ops.equalDOF(3, 4, 2)       # top face stays horizontal (node 3 retained)
    # nDMaterial PM4Sand tag Dr G0 hpo rho  P_atm h0 emax emin nb nd Ado z_max cz ce phi_cv nu
    ops.nDMaterial("PM4Sand", 1, pr["Dr"], pr["G0"], pr["hpo"], 0.0,
                   pr["p_a"], -1.0, E_MAX, E_MIN, 0.5, 0.1, -1.0, -1.0, 250.0, -1.0, pr["phi_cv"], pr["nu"])
    ops.element("SSPquad", 1, 1, 2, 3, 4, 1, "PlaneStrain", 1.0, 0.0, 0.0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, -p0, -p0)       # retained corner carries each unit face's force: sigma_xx = sigma_yy = -p0
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    # 1e-4 p0 (0.01 kPa at 100 kPa): PM4Sand integrates explicitly inside the step, and past ~2 % axial strain
    # Newton stalls at a residual of ~1e-3 kPa with the displacement correction already ~1e-8; K4 compares at
    # the percent level, so a 1e-4 relative force balance is ample.
    ops.test("NormUnbalance", 1.0e-4 * p0, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")
    if ops.analyze(20) != 0:
        return dict(status="consolidation_failed")
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.setParameter("-val", 0, "-ele", 1, "FirstCall", "1")   # re-initialise PM4Sand at the consolidated stress
    ops.integrator("LoadControl", 0.0)
    if ops.analyze(1) != 0:     # one step at constant load (factor stays 1): stage switch settles
        return dict(status="stage_switch_failed")
    ops.loadConst("-time", 0.0)

    ux0, uy0 = ops.nodeDisp(3, 1), ops.nodeDisp(3, 2)
    d = eps_end / n
    if kind == "PS":
        ops.timeSeries("Linear", 2)
        ops.pattern("Plain", 2, 2)
        ops.load(3, 0.0, -1.0)  # DisplacementControl's reference load (the lateral force stays in the constant pattern)

        def set_step(frac):
            ops.integrator("DisplacementControl", 3, 2, -d * frac)
    elif kind == "PSU":
        # -useLast: a Path series returns 0 past its last time, and n steps of 1/n can overshoot 1.0 by round-off
        ops.timeSeries("Path", 2, "-time", 0.0, 1.0, "-values", uy0, uy0 - eps_end, "-useLast")
        ops.timeSeries("Path", 3, "-time", 0.0, 1.0, "-values", ux0, ux0 + eps_end, "-useLast")
        ops.pattern("Plain", 2, 2)
        ops.sp(3, 2, 1.0)
        ops.pattern("Plain", 3, 3)
        ops.sp(3, 1, 1.0)

        def set_step(frac):
            ops.integrator("LoadControl", frac / n)
    else:
        return dict(status=f"kind {kind} not plane strain: PM4Sand is a 2D model")
    set_step(1.0)
    n_sub = [0]

    def advance():
        """One output step. Fallbacks, in order: Newton; KrylovNewton; NewtonLineSearch; then 10 substeps with
        KrylovNewton. Counts the steps that needed substepping (reported, not hidden)."""
        for alg in (("Newton",), ("KrylovNewton",), ("NewtonLineSearch", 0.8)):
            ops.algorithm(*alg)
            if ops.analyze(1) == 0:
                ops.algorithm("Newton")
                return True
        ops.algorithm("KrylovNewton")
        set_step(0.1)
        ok = all(ops.analyze(1) == 0 for _ in range(10))
        set_step(1.0)
        ops.algorithm("Newton")
        n_sub[0] += 1
        return ok

    rows, status = [], "ok"

    def record():
        sx, sy, txy = ops.eleResponse(1, "stress")[:3]
        exx = ops.nodeDisp(3, 1)            # unit square: displacement = strain
        eyy = ops.nodeDisp(3, 2)
        rows.append(dict(exx=exx, eyy=eyy, sxx=sx, syy=sy, txy=txy))

    record()
    for k in range(n):
        if not advance():
            status = f"diverged at step {k + 1}"
            break
        record()
        if not all(math.isfinite(v) for v in rows[-1].values()):
            rows.pop()
            status = f"nan at step {k + 1}"   # analyze() returned 0 with a NaN stress: stop, do not report it
            break
    # strains reported from the start of shearing (the consolidation strain is not part of the test)
    r0 = dict(rows[0])
    for r in rows:
        r["exx"] -= r0["exx"]
        r["eyy"] -= r0["eyy"]
    return dict(status=status, params=pr, rows=rows, substepped_steps=n_sub[0])


def main(argv=None):
    """Each case runs in its own child process: two PM4Sand cases in one interpreter interfere even across wipe()
    (the second consolidation returns NaN; each case alone is clean), so no state may carry between cases."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--dist", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--one", default=None, help="internal: run only this case id and print its JSON")
    a = ap.parse_args(argv)
    spec = json.load(open(a.cases, encoding="utf-8"))
    if a.one is not None:
        os.add_dll_directory(a.dist)
        sys.path.insert(0, a.dist)
        import opensees as ops
        case = next(c for c in spec["cases"] if c["id"] == a.one)
        res = run_case(ops, case, spec["sand"], a.n)
        res["build"] = ops.ladrunoBuild()
        sys.stdout.write("K4JSON " + json.dumps(res) + "\n")
        return
    import subprocess
    out = dict(sand=spec["sand"], results={})
    for case in spec["cases"]:
        if case["kind"] not in ("PS", "PSU"):
            continue
        p = subprocess.run([sys.executable, "-S", os.path.abspath(__file__), "--dist", a.dist, "--cases", a.cases,
                            "--out", a.out, "--n", str(a.n), "--one", case["id"]],
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        line = next((l for l in p.stdout.splitlines() if l.startswith("K4JSON ")), None)
        res = json.loads(line[7:]) if line else dict(status=f"child failed rc={p.returncode}: {p.stderr[-300:]}")
        out["build"] = res.pop("build", out.get("build"))
        out["results"][case["id"]] = res
        print(case["id"], res["status"], flush=True)
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(out, f)


if __name__ == "__main__":
    main()
