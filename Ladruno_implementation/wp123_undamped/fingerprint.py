#!/usr/bin/env python3
"""WP-123 bit-identity fingerprint: record every value the four coupling elements produce
in damped dynamics, as exact float reprs, so two builds can be compared byte for byte.

Scenarios (per element): Newmark and HHT with rayleigh ON and OFF, 40 steps, disp + vel
of every free DOF each step, then dampingForce / force responses; explicit
CentralDifferenceLadruno with and without betaK: self-reported dtcr, integrator
criticalTimeStep, 40 steps of disp. Uses the models of tests/test_ladruno_undamped_couplings.py.

    <py3.12> -S fingerprint.py OUT.json        (same -S bootstrap as run_pytest.py)
    <py3.12> -S fingerprint.py --compare A.json B.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))


def compare(a, b):
    A = json.load(open(a, encoding="utf-8"))
    B = json.load(open(b, encoding="utf-8"))
    keys = sorted(set(A) | set(B))
    diff = [k for k in keys if A.get(k) != B.get(k)]
    nvals = sum(len(v) for v in A.values())
    print(f"{len(keys)} series, {nvals} recorded values; {len(diff)} series differ")
    for k in diff[:20]:
        print("  DIFF", k)
    return 1 if diff else 0


def main():
    if sys.argv[1] == "--compare":
        return compare(sys.argv[2], sys.argv[3])
    DIST = os.path.join(ROOT, "dist", "bin")
    assert sys.flags.no_site, "run with python -S"
    os.add_dll_directory(DIST)
    sys.path.insert(0, DIST)
    sys.path.insert(0, os.path.join(ROOT, "tests"))
    sys.path.append(os.path.join(os.path.dirname(sys.executable), "Lib", "site-packages"))  # pytest
    os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
    import opensees
    assert os.path.normcase(opensees.__file__) == os.path.normcase(os.path.join(DIST, "opensees.pyd"))
    import test_ladruno_undamped_couplings as T
    ops = T.ops
    out = {}
    for name in T.ELEMENTS:
        for integ in sorted(T.INTEGRATORS):
            for label, ray in (("on", T.RAYLEIGH_ON), ("off", T.RAYLEIGH_OFF)):
                hist, damp = T._implicit_run(name, ray, integ, nsteps=40)
                out[f"{name}/{integ}/{label}/traj"] = [repr(x) for step in hist for dv in step for x in dv]
                out[f"{name}/{integ}/{label}/dampingForce"] = [repr(x) for x in damp]
                out[f"{name}/{integ}/{label}/force"] = [repr(x) for x in ops.eleResponse(1, "force")]
        for label, ray in (("betaK", (0.0, 1.0e-3, 0.0, 0.0)), ("off", T.RAYLEIGH_OFF)):
            ops.wipe()
            free = T.MODELS[name]()
            ops.rayleigh(*ray)
            rec = [repr(ops.eleResponse(1, "dtcr")[0])]
            ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("Diagonal")
            ops.test("NormDispIncr", 1.0e-12, 1); ops.algorithm("Linear")
            ops.integrator("CentralDifferenceLadruno", "-cfl"); ops.analysis("Transient")
            for _ in range(40):
                assert ops.analyze(1, 1.0e-5) == 0
                rec += [repr(ops.nodeDisp(n, d)) for n, d in free]
            rec.append(repr(ops.criticalTimeStep()))
            out[f"{name}/explicit/{label}"] = rec
    json.dump(out, open(sys.argv[1], "w", encoding="utf-8"), indent=0)
    print(f"wrote {len(out)} series, {sum(len(v) for v in out.values())} values -> {sys.argv[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
