"""WP-240: before/after evidence for the single element update per DR step.

Run it once with each build (a pre-WP-240 `dist` and the WP-240 one) and
compare the two JSON files:

    py -3.12 dr_single_update_bench.py run <dist/bin> <out.json>
    py -3.12 dr_single_update_bench.py compare <before.json> <after.json>

`run` records
  * the observables of tests/test_wp240_dr_single_update.py (the `_REFERENCE`
    table of that file is the `reference` block of the BEFORE run);
  * for three larger models, every nodal displacement after the march as hex
    floats (bit-exact comparison) and the wall time per step, best of 3:
      brick_j2      16x16x16 LadrunoBrick, LadrunoJ2, top face driven by `sp`
                    under the Transformation handler (32768 Gauss points)
      stdbrick_j2   16x16x8 stdBrick, J2Plasticity, same loading (16384 GP)
      quad_j2       64x64 LadrunoQuad plane strain, LadrunoJ2, edge load
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))

E, NU, SIG0 = 2.0e5, 0.3, 250.0
K, G = E / (3.0 * (1 - 2 * NU)), E / (2.0 * (1 + NU))


def _load_engine(dist_bin):
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(dist_bin)
    os.environ["PATH"] = dist_bin + os.pathsep + os.environ.get("PATH", "")
    sys.path.insert(0, dist_bin)
    import opensees as ops
    got = os.path.dirname(os.path.abspath(ops.__file__))
    assert os.path.normcase(got) == os.path.normcase(dist_bin), ops.__file__
    return ops


def _brick(ops, nx, ny, nz, element, rate):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    if element == "LadrunoBrick":
        ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", SIG0, 0.0, 0.0, 200.0)
    else:
        ops.nDMaterial("J2Plasticity", 1, K, G, SIG0, SIG0, 0.0, 200.0)
    h = 1.0 / nx
    tag = {}
    t = 1
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(t, i * h, j * h, k * h)
                tag[(i, j, k)] = t
                t += 1
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                c = [tag[(i, j, k)], tag[(i + 1, j, k)], tag[(i + 1, j + 1, k)],
                     tag[(i, j + 1, k)], tag[(i, j, k + 1)], tag[(i + 1, j, k + 1)],
                     tag[(i + 1, j + 1, k + 1)], tag[(i, j + 1, k + 1)]]
                ops.element(element, e, *c, 1)
                e += 1
    for (i, j, k), nd in tag.items():
        if k == 0:
            ops.fix(nd, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for (i, j, k), nd in tag.items():
        if k == nz:
            ops.sp(nd, 3, -1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1.0e-12, 1, 0)
    ops.algorithm("Linear")
    ops.integrator("LadrunoDynamicRelaxation", "-dt", rate)
    ops.analysis("Transient")
    return sorted(tag.values())


def _quad(ops, n):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", SIG0, 0.0, 0.0, 200.0)
    h = 1.0 / n
    tag = {}
    t = 1
    for j in range(n + 1):
        for i in range(n + 1):
            ops.node(t, i * h, j * h)
            tag[(i, j)] = t
            t += 1
    e = 1
    for j in range(n):
        for i in range(n):
            ops.element("LadrunoQuad", e, tag[(i, j)], tag[(i + 1, j)],
                        tag[(i + 1, j + 1)], tag[(i, j + 1)], 1,
                        "-type", "PlaneStrain", "-thick", 1.0)
            e += 1
    for j in range(n + 1):
        ops.fix(tag[(0, j)], 1, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(n + 1):
        w = 0.5 if j in (0, n) else 1.0
        ops.load(tag[(n, j)], 1.2 * SIG0 * w / n, 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1.0e-12, 1, 0)
    ops.algorithm("Linear")
    ops.integrator("LadrunoDynamicRelaxation")
    ops.analysis("Transient")
    return sorted(tag.values())


_CASES = {
    # name: (builder, steps, dt). The sp rate reaches 3x the yield strain.
    "brick_j2": (lambda ops: _brick(ops, 16, 16, 16, "LadrunoBrick", 1.0e-5), 400, 1.0e-5),
    "stdbrick_j2": (lambda ops: _brick(ops, 16, 16, 8, "stdBrick", 1.0e-5), 400, 1.0e-5),
    "quad_j2": (lambda ops: _quad(ops, 64), 400, 1.0),
}


def _march(ops, name):
    build, steps, dt = _CASES[name]
    best = None
    state = None
    for _ in range(3):
        nodes = build(ops)
        assert ops.analyze(1, dt) == 0          # domainChanged + the M* build
        t0 = time.perf_counter()
        assert ops.analyze(steps - 1, dt) == 0
        el = (time.perf_counter() - t0) / (steps - 1)
        best = el if best is None else min(best, el)
        got = [float(u).hex() for nd in nodes for u in ops.nodeDisp(nd)]
        got.append(float(ops.ladrunoDR("residualNorm")).hex())
        got.append(float(ops.ladrunoDR("kineticEnergy")).hex())
        assert state is None or state == got, name + ": not reproducible run to run"
        state = got
    ops.wipe()
    return {"steps": steps, "ms_per_step": 1.0e3 * best, "state": state}


def run(dist_bin, out):
    ops = _load_engine(os.path.abspath(dist_bin))
    sys.path.insert(0, os.path.join(ROOT, "tests"))
    import test_wp240_dr_single_update as t
    t.ops = ops                                  # the engine under test, whatever _testbed found
    ref = {
        "brick_kinetic": t._run_brick(t._KINETIC),
        "brick_viscous": t._run_brick(t._VISCOUS),
        "quad_kinetic": t._run_quad(t._KINETIC),
        "quad_viscous": t._run_quad(t._VISCOUS),
    }
    build = None
    try:
        build = str(ops.ladrunoBuild())
    except Exception:
        pass
    res = {"build": build, "reference": ref, "cases": {}}
    for name in _CASES:
        res["cases"][name] = _march(ops, name)
        print("%-12s %8.3f ms/step" % (name, res["cases"][name]["ms_per_step"]), flush=True)
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)


def compare(before, after):
    with open(before) as fh:
        a = json.load(fh)
    with open(after) as fh:
        b = json.load(fh)
    print("before: %s" % (a.get("build"),))
    print("after : %s" % (b.get("build"),))
    ok = True
    for name, ref in a["reference"].items():
        for key, v in ref.items():
            w = b["reference"][name][key]
            same = float(v).hex() == float(w).hex()
            ok = ok and same
            print("  %-14s %-10s %-24.17g %-24.17g %s" % (
                name, key, v, w, "bit-identical" if same else "DIFFERENT"))
    for name, ca in a["cases"].items():
        cb = b["cases"][name]
        n = len(ca["state"])
        diff = sum(1 for x, y in zip(ca["state"], cb["state"]) if x != y)
        ok = ok and diff == 0 and n == len(cb["state"])
        print("%-12s %6d values, %d differ | %8.3f -> %8.3f ms/step (x%.2f)" % (
            name, n, diff, ca["ms_per_step"], cb["ms_per_step"],
            ca["ms_per_step"] / cb["ms_per_step"]))
    print("SAME ANSWER" if ok else "ANSWERS DIFFER")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "compare":
        sys.exit(compare(sys.argv[2], sys.argv[3]))
    else:
        sys.exit(__doc__)
