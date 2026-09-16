"""WP-107 (ADR-75b L3-1) -- the threaded-Domain::update bench and identity gate.

A plane-strain LadrunoQuad strip footing, the TIMs-Workbench shape: LadrunoQuad
`-formulation bbar`, half a strip footing on a SANISAND (CP1/Gorini) half-space,
self-weight, prescribed settlement, global Newton. Derived from
`../../hypo_bearing/adr92_f7/strip_quad_implex.py`, with three changes that matter:

  * `-implex` is OFF. WP-107 REFUSES to thread a LadrunoSANISAND running under
    -implex, because its diagnostics are a process-wide ledger with FLOATING-POINT
    accumulators (see Ladruno_implementation/107_ladruno_openmp_element_loop.md
    section 5). The `--implex` switch here exists to demonstrate that refusal.
  * the mesh is refinable (`--h`), because the whole question is what happens at
    tens of thousands of Gauss points;
  * it reports wall time per step and writes the load-settlement curve at FULL
    double precision (`repr`), because the acceptance gate is BIT-IDENTITY across
    thread counts, not a plot that looks the same.

Also carries an elastic arm (`--mat elastic`) whose material,
ElasticIsotropicPlaneStrain2D, is separately allowlisted -- a non-SANISAND
bit-identity check on the same element.

Usage (from the worktree root):

    set PYTHONPATH=<worktree>\\dist\\bin
    python3.12 -u Ladruno_files/testbed/perf/wp107/wp107_strip_bench.py \\
        --threads 4 --h 0.15 --steps 12 --out out_t4.csv

Acceptance protocol (ADR-75b section 7, correctness protocol item 4): run each
thread count several times and require EVERY run to reproduce itself AND the
serial baseline, with MKL_NUM_THREADS pinned to 1 so the solver's own ~1 ULP
jitter (ADR-75b section 3 P-6) cannot masquerade as an assembly race.
"""
import argparse
import os
import sys
import time

# MKL_NUM_THREADS must be pinned BEFORE the module loads, or the solver's own
# run-to-run jitter makes the bit-identity gate unfalsifiable (ADR-75b P-6).
os.environ.setdefault("MKL_NUM_THREADS", "1")

import opensees as ops   # noqa: E402  (PYTHONPATH must point at a dist/bin)

# --- CP1 / Gorini parameter set -- tests/test_ladruno_sanisand.py::_PARAMS ----
_PARAMS = [
    264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0,
    0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0,
]
_RHO = _PARAMS[17]          # 2.0 Mg/m3 with kPa/m/s units
_G = 9.81

_B_HALF = 1.0               # half footing width [m]
_XMAX = 6.0
_ZMIN = -6.0

# Substep cap 1000, as the deck this is derived from uses. UNCAPPED (0) is the
# vanilla default and it HANGS this deck: a badly conditioned integration point
# substeps down to dT_min = 1e-6, i.e. up to a million substeps at ONE point.
_MAXSUBSTEPS = 1000
_PMIN = 0.0101

_DS = 0.0005                # prescribed settlement increment [m/step]
_TOL = 1.0e-4
_MAXITER = 40
_N_GRAV = 5

_MAT_TAG = 9902


def _grid(h):
    nx = int(round(_XMAX / h))
    nz = int(round(-_ZMIN / h))
    nodes = {}
    tag = 1
    for i in range(nx + 1):
        for k in range(nz + 1):
            nodes[(i, k)] = tag
            ops.node(tag, i * h, _ZMIN + k * h)
            tag += 1
    return nx, nz, nodes


def build(h, mat, scheme, tangent, implex):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    nx, nz, nodes = _grid(h)

    if mat == "sanisand":
        # IntScheme / TanType / JacoType / TolF / TolR are POSITIONAL in this
        # parser (they follow Rho), not flags. The flag options come after.
        pos = [scheme, tangent, 1, 1.0e-7, 1.0e-7]
        opts = ["-Pmin", _PMIN, "-maxSubsteps", _MAXSUBSTEPS]
        if implex:
            opts += ["-implex"]
        ops.nDMaterial("LadrunoSANISAND", _MAT_TAG, *_PARAMS, *pos, *opts)
    elif mat == "elastic":
        ops.nDMaterial("ElasticIsotropic", _MAT_TAG, 60000.0, 0.3, _RHO)
    else:
        raise SystemExit("unknown --mat %s" % mat)

    eid = 1
    for i in range(nx):
        for k in range(nz):
            ops.element("LadrunoQuad", eid,
                        nodes[(i, k)], nodes[(i + 1, k)],
                        nodes[(i + 1, k + 1)], nodes[(i, k + 1)],
                        _MAT_TAG,
                        "-thick", 1.0, "-type", "PlaneStrain",
                        "-formulation", "bbar",
                        "-rho", _RHO, "-body", 0.0, -_G)
            eid += 1

    mask = {}
    for k in range(nz + 1):
        mask[nodes[(0, k)]] = [1, 0]
        mask[nodes[(nx, k)]] = [1, 0]
    for i in range(nx + 1):
        m = mask.setdefault(nodes[(i, 0)], [0, 0])
        m[1] = 1
    for n, m in mask.items():
        ops.fix(n, m[0], m[1])

    footing = [nodes[(i, nz)] for i in range(nx + 1) if i * h <= _B_HALF + 1e-9]
    top = [nodes[(i, nz)] for i in range(nx + 1)]
    return nodes, footing, top, eid - 1


def _analysis(dlambda, system):
    # The SANISAND tangent is UNSYMMETRIC (non-associated flow).
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system(system)
    ops.test("NormDispIncr", _TOL, _MAXITER, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", dlambda)
    ops.analysis("Static")


def main(argv=None):
    global _DS, _MAXSUBSTEPS, _TOL
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--h", type=float, default=0.25, help="element size [m]")
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--mat", default="sanisand", choices=("sanisand", "elastic"))
    ap.add_argument("--scheme", type=int, default=1, help="SANISAND IntScheme")
    ap.add_argument("--tangent", type=int, default=0, help="SANISAND TanType")
    ap.add_argument("--implex", action="store_true",
                    help="turn -implex ON (WP-107 then REFUSES to thread; that is the point)")
    ap.add_argument("--system", default="BandGeneral")
    ap.add_argument("--load", default="oedometer",
                    choices=("oedometer", "footing"))
    ap.add_argument("--profile", default=None,
                    help="write a DEEP profiler report here. Deep profiling and "
                         "threading are mutually exclusive by design (hazard H6), "
                         "so use this on the 1-thread baseline only.")
    ap.add_argument("--ds", type=float, default=_DS)
    ap.add_argument("--maxsub", type=int, default=_MAXSUBSTEPS)
    ap.add_argument("--tol", type=float, default=_TOL)
    ap.add_argument("--out", default="wp107_bench.csv")
    args = ap.parse_args(argv)

    _DS, _MAXSUBSTEPS, _TOL = args.ds, args.maxsub, args.tol

    print("build:", ops.ladrunoBuild().strip().splitlines()[0])
    got = ops.ladrunoThreads(args.threads)
    print("ladrunoThreads requested=%d stored=%s  MKL_NUM_THREADS=%s"
          % (args.threads, got, os.environ.get("MKL_NUM_THREADS")))

    nodes, footing, top, nele = build(args.h, args.mat, args.scheme,
                                      args.tangent, args.implex)
    # WHY THE DEFAULT LOAD CASE IS THE CONFINED ONE. A strip FOOTING on a free
    # surface drives the corner Gauss points into tension (p < 0), which is the
    # zero-confinement wall ADR-93 is about: the material clamps to p_min,
    # ModifiedEuler substeps to its cap, and the step refuses. That is a real and
    # interesting failure -- it is just not an instrument. This WP has to measure
    # element-loop THROUGHPUT and BIT-IDENTITY on the same element + material, so
    # the default is a laterally confined compression (`oedometer`): rollers on
    # both sides, pinned base, uniform prescribed settlement on the WHOLE top
    # surface. p rises monotonically, every Gauss point integrates, nothing
    # refuses, and the kernel exercised is identical. `--load footing` keeps the
    # original path for anyone who wants it.
    driven = top if args.load == "oedometer" else footing
    ngp = nele * 4
    print("elements: %d   Gauss points: %d   driven nodes: %d   dof ~ %d   load=%s"
          % (nele, ngp, len(driven), 2 * len(nodes), args.load))

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    _analysis(1.0 / _N_GRAV, args.system)
    if args.mat == "sanisand":
        ops.updateMaterialStage("-material", _MAT_TAG, "-stage", 0)
    t_grav = time.perf_counter()
    for s in range(_N_GRAV):
        rc = ops.analyze(1)
        if rc != 0:
            raise SystemExit("gravity step %d failed (rc=%d)" % (s + 1, rc))
    t_grav = time.perf_counter() - t_grav
    ops.loadConst("-time", 0.0)
    if args.mat == "sanisand":
        ops.updateMaterialStage("-material", _MAT_TAG, "-stage", 1)

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for n in driven:
        ops.sp(n, 2, -1.0)
    ops.integrator("LoadControl", _DS)

    if args.profile:
        ops.profiler("start", "-deep", "-perStep")

    rows = []
    t_steps = []
    for s in range(args.steps):
        t0 = time.perf_counter()
        rc = ops.analyze(1)
        t_steps.append(time.perf_counter() - t0)
        ops.reactions()
        q = sum(ops.nodeReaction(n, 2) for n in driven)
        w = -ops.nodeDisp(driven[0], 2)
        rows.append((s + 1, rc, w, -q))
        if rc != 0:
            print("STOPPED at step %d (rc=%d)" % (s + 1, rc))
            break

    if args.profile:
        ops.profiler("stop")
        ops.profiler("report", args.profile)
        print("profiler report ->", args.profile)

    # repr(), not %g: the gate is BIT-identity, so every bit has to survive the
    # round trip through the file.
    with open(args.out, "w", newline="\n") as fh:
        fh.write("step,rc,settlement_m,footing_load_kN_per_m\n")
        for r in rows:
            fh.write("%d,%d,%s,%s\n" % (r[0], r[1], repr(r[2]), repr(r[3])))

    tot = sum(t_steps)
    print("gravity wall  : %.3f s (%d steps)" % (t_grav, _N_GRAV))
    print("settle wall   : %.3f s over %d steps" % (tot, len(t_steps)))
    print("per-step wall : %.4f s (mean), %.4f s (min), %.4f s (median)"
          % (tot / max(1, len(t_steps)), min(t_steps) if t_steps else 0.0,
             sorted(t_steps)[len(t_steps) // 2] if t_steps else 0.0))
    print("per-step list : " + " ".join("%.4f" % t for t in t_steps))
    print("wrote %s (%d rows)" % (args.out, len(rows)))
    print("WALLLINE threads=%d h=%s mat=%s scheme=%d tan=%d nele=%d "
          "grav=%.4f settle=%.4f perstep=%.5f median=%.5f"
          % (args.threads, args.h, args.mat, args.scheme, args.tangent, nele,
             t_grav, tot, tot / max(1, len(t_steps)),
             sorted(t_steps)[len(t_steps) // 2] if t_steps else 0.0))


if __name__ == "__main__":
    main()
