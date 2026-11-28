"""M1 harness (WP-263, decision D2): does the exact RC tangent cut global Newton iterations?

Runs existing LadrunoRCConcrete decks twice -- the default tangent (mode 0) and `-dualTangent`
(mode 3, one Dual<6> pass of the return map) -- with the SAME solver settings, and prints, per
case and tangent: Newton iterations per step, total iterations, failed steps (and the fallback
that rescued them), wall time, and `ladrunoBuild()`.

Decks (geometry, backbones and steel are tests/_testbed/rc_wall_harness.py's, imported, not
copied):
  panel_implicit  a 2x2 ASDShellQ4 / LayeredShell RC panel (4 concrete + 2 smeared Steel02 web
                  layers), base fixed, top edge tied to one control node, monotonic in-plane shear
                  by DisplacementControl to 1.5 mm. Implicit material: the consistent-vs-secant
                  question in its plain form.
  panel_implex    the same panel with -implex (the wall harness's robustness mode).
  wall_implex     the 4x4 squat-wall deck of rc_wall_harness.py (-beta -lublinerReduced -interlock
                  -cyclic -implex), its first cycles: amplitudes 0.3, 0.3, 0.6 mm.
The wall stalls past first crack for a STRUCTURAL reason (quirk
LEGACY-cyclic-softening-rc-shell-wall-asdshellq4.md): a stall there is a finding to report, not
something the tangent is expected to fix.

Solver: `test NormDispIncr <tol> <maxIter>` (default 1e-8 mm, 200), `algorithm Newton` (full
Newton, the tangent re-formed every iteration -- what measures tangent quality; pass
--algorithm KrylovNewton for the wall harness's own setting). A step that fails is retried with
the wall harness's fallbacks (NewtonLineSearch, ModifiedNewton, Broyden) and counted as failed:
its entry in iters/step is the failed attempt's count (maxIter + 1), the rescue's own count is
listed under "rescued"; a step no fallback rescues ends that run.

Run after building, against the freshly built engine (the orchestrator; WP-263 does not build):

    py -3.12 -S tests/_testbed/rc_newton_count.py --dist dist/bin [--out m1.json]
    py -3.12 -S tests/_testbed/rc_newton_count.py --dry-run          # no engine needed

`-S` keeps a site-packages .pth from importing another build's opensees (WP-175/176); the harness
needs only the standard library and the engine.
"""
import argparse
import json
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_TESTS = os.path.dirname(_HERE)
_REPO = os.path.dirname(_TESTS)

ALGOS_FALLBACK = ("NewtonLineSearch", "ModifiedNewton", "Broyden")


def _wall_constants():
    """The deck constants of rc_wall_harness.py, read without importing an engine."""
    import ast
    import operator
    ops_ = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}

    def ev(n, env):   # numeric literals, + - * /, unary -, lists/tuples, names already read
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)):
            return n.value
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return -ev(n.operand, env)
        if isinstance(n, ast.BinOp) and type(n.op) in ops_:
            return ops_[type(n.op)](ev(n.left, env), ev(n.right, env))
        if isinstance(n, (ast.List, ast.Tuple)):
            return [ev(e, env) for e in n.elts]
        if isinstance(n, ast.Name) and n.id in env:
            return env[n.id]
        raise ValueError(ast.dump(n))

    src = open(os.path.join(_HERE, "rc_wall_harness.py"), encoding="utf-8").read()
    out = {}
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign):
            try:
                val = ev(node.value, out)
            except ValueError:
                continue
            targets = node.targets[0].elts if isinstance(node.targets[0], ast.Tuple) else node.targets
            vals = val if isinstance(node.targets[0], ast.Tuple) else (val,)
            for t, v in zip(targets, vals):
                if isinstance(t, ast.Name):
                    out[t.id] = v
    return out


W = _wall_constants()

CASES = {
    "panel_implicit": dict(nx=2, ny=2, flags=("-beta", "-lublinerReduced"), implex=False,
                           protocol=("mono", 1.5, 60)),
    "panel_implex":   dict(nx=2, ny=2, flags=("-beta", "-lublinerReduced"), implex=True,
                           protocol=("mono", 1.5, 60)),
    "wall_implex":    dict(nx=W["NX"], ny=W["NY"],
                           flags=("-beta", "-lublinerReduced", "-interlock", "-cyclic",
                                  "-crackSpacing", 80.0, "-agg", 16.0),
                           implex=True, protocol=("cyclic", (0.3, 0.3, 0.6), 80)),
}
MODES = ("default", "dual")


def material_args(case, mode):
    c = CASES[case]
    a = ["LadrunoRCConcrete", 1, W["E"], W["NU"],
         "-Ce", *W["CE"], "-Cs", *W["CS"], "-Cd", *W["CD"],
         "-Te", *W["TE"], "-Ts", *W["TS"], "-Td", *W["TD"], *c["flags"]]
    if c["implex"]:
        a.append("-implex")
    if mode == "dual":
        a.append("-dualTangent")
    return a


def targets(case):
    kind, a, n = CASES[case]["protocol"]
    if kind == "mono":
        return [(a, n)]
    return [(t, n) for amp in a for t in (amp, -amp)]


def build(ops, case, mode):
    c = CASES[case]
    nx, ny = c["nx"], c["ny"]
    L, H, TH = W["L"], W["H"], W["TH"]
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nid, k = {}, 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            nid[(i, j)] = k
            ops.node(k, i * L / nx, j * H / ny, 0.0)
            k += 1
    ops.nDMaterial(*material_args(case, mode))
    ops.uniaxialMaterial("Steel02", 11, W["Fy"], W["Es"], W["BSTEEL"], 18.0, 0.925, 0.15)
    ops.nDMaterial("PlateRebar", 21, 11, 90.0)
    ops.nDMaterial("PlateRebar", 22, 11, 0.0)
    tc = TH * 0.94
    ops.section("LayeredShell", 1, 6, 1, tc / 4, 1, tc / 4, 1, tc / 4, 1, tc / 4,
                21, W["RHO_V"] * TH, 22, W["RHO_H"] * TH)
    eid = 1
    for j in range(ny):
        for i in range(nx):
            ops.element("ASDShellQ4", eid, nid[(i, j)], nid[(i + 1, j)], nid[(i + 1, j + 1)], nid[(i, j + 1)], 1)
            eid += 1
    base = [nid[(i, 0)] for i in range(nx + 1)]
    top = [nid[(i, ny)] for i in range(nx + 1)]
    for n in base:
        ops.fix(n, 1, 1, 1, 1, 1, 1)
    for (i, j), n in nid.items():
        if j > 0:
            ops.fix(n, 0, 0, 1, 1, 1, 0)
    ctrl = top[nx // 2]
    for n in top:
        if n != ctrl:
            ops.equalDOF(ctrl, n, 1)
    return ctrl, base


def run_one(ops, case, mode, tol, maxit, algorithm):
    ctrl, base = build(ops, case, mode)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(ctrl, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", tol, maxit, 0)
    ops.algorithm(algorithm)
    ops.analysis("Static")
    iters, failed, rescued, curve = [], [], [], []
    cur, step, stopped = 0.0, 0, None
    t0 = time.perf_counter()
    for target, n in targets(case):
        ops.integrator("DisplacementControl", ctrl, 1, (target - cur) / n)
        for _ in range(n):
            step += 1
            ok = ops.analyze(1)
            it = ops.testIter()
            if ok != 0:
                failed.append(step)
                for alg in ALGOS_FALLBACK:
                    ops.algorithm(alg)
                    ok = ops.analyze(1)
                    ops.algorithm(algorithm)
                    if ok == 0:
                        rescued.append((step, alg, ops.testIter()))
                        break
                if ok != 0:
                    stopped = (step, target, cur)
                    break
            iters.append(it)
            ops.reactions()
            curve.append((ops.nodeDisp(ctrl, 1), -sum(ops.nodeReaction(b, 1) for b in base)))
        if stopped:
            break
        cur = target
    wall = time.perf_counter() - t0
    return dict(case=case, mode=mode, steps=len(iters), iters=iters, total=sum(iters),
                mean=(sum(iters) / len(iters)) if iters else float("nan"),
                max=max(iters) if iters else 0, failed=failed, rescued=rescued,
                stopped=stopped, wall_s=wall, curve=curve)


def import_engine(dist):
    if dist:
        dist = os.path.abspath(dist)
        if hasattr(os, "add_dll_directory"):
            os.add_dll_directory(dist)
        sys.path.insert(0, dist)
    import opensees as ops   # the engine under test (no openseespy fallback: a wrong engine must fail)
    return ops


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dist", default=None, help="directory holding opensees.pyd (e.g. dist/bin)")
    ap.add_argument("--cases", default=",".join(CASES), help="comma list of " + ", ".join(CASES))
    ap.add_argument("--tol", type=float, default=1.0e-8)
    ap.add_argument("--maxit", type=int, default=200)
    ap.add_argument("--algorithm", default="Newton")
    ap.add_argument("--modes", default=",".join(MODES), help="comma list of default, dual")
    ap.add_argument("--out", default=None, help="write the full result (per-step iterations) as JSON")
    ap.add_argument("--dry-run", action="store_true", help="print the plan; import no engine")
    args = ap.parse_args(argv)
    cases = [c.strip() for c in args.cases.split(",") if c.strip()]
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    for m in modes:
        if m not in MODES:
            ap.error(f"unknown mode {m!r}")
    for c in cases:
        if c not in CASES:
            ap.error(f"unknown case {c!r}")

    if args.dry_run:
        print(f"rc_newton_count dry run: test NormDispIncr {args.tol:g} {args.maxit}, algorithm {args.algorithm}")
        for c in cases:
            steps = sum(n for _, n in targets(c))
            print(f"  {c}: {CASES[c]['nx']}x{CASES[c]['ny']} ASDShellQ4, {steps} steps, targets "
                  f"{[t for t, _ in targets(c)]} mm")
            for m in modes:
                print(f"    {m:7s}: nDMaterial " + " ".join(str(x) for x in material_args(c, m)))
        return 0

    ops = import_engine(args.dist)
    build_sha = ops.ladrunoBuild() if hasattr(ops, "ladrunoBuild") else "unknown"
    print(f"rc_newton_count: engine {getattr(ops, '__file__', '?')}  ladrunoBuild {build_sha}")
    print(f"solver: test NormDispIncr {args.tol:g} {args.maxit}, algorithm {args.algorithm}")
    results = []
    for c in cases:
        for m in modes:
            r = run_one(ops, c, m, args.tol, args.maxit, args.algorithm)
            results.append(r)
            print(f"{c:15s} {m:7s} steps {r['steps']:4d}  iters {r['total']:6d}  mean {r['mean']:6.2f}  "
                  f"max {r['max']:4d}  failed {len(r['failed']):3d}  wall {r['wall_s']:8.2f} s"
                  + (f"  STOPPED at step {r['stopped'][0]} (target {r['stopped'][1]} mm)" if r["stopped"] else ""))
            print(f"    iters/step: {' '.join(str(i) for i in r['iters'])}")
            if r["rescued"]:
                print(f"    rescued: {r['rescued']}")
    print("\nsummary (dual / default):")
    for c in cases:
        d0 = next((r for r in results if r["case"] == c and r["mode"] == "default"), None)
        d3 = next((r for r in results if r["case"] == c and r["mode"] == "dual"), None)
        if d0 is None or d3 is None:
            continue
        n = min(d0["steps"], d3["steps"])
        same0, same3 = sum(d0["iters"][:n]), sum(d3["iters"][:n])
        print(f"  {c:15s} iterations over the {n} common steps: {same3}/{same0} = "
              f"{(same3 / same0) if same0 else float('nan'):.3f}; wall {d3['wall_s']:.2f}/{d0['wall_s']:.2f} s; "
              f"failed {len(d3['failed'])}/{len(d0['failed'])}; steps {d3['steps']}/{d0['steps']}")
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(dict(build=build_sha, tol=args.tol, maxit=args.maxit, algorithm=args.algorithm,
                           date=time.strftime("%Y-%m-%d"), results=results), fh, indent=1)
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
