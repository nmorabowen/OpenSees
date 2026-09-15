"""WP-99 (F7) -- the plane-strain LadrunoQuad strip that reproduces the silent
commit-time companion failure, and its cure.

A deliberately COARSE plane-strain twin of the ADR-95 hypo_bearing slab deck:
half a strip footing on a SANISAND (CP1/Gorini) half-space, self-weight on,
`-implex` WITHOUT `-implexControl`, `-maxSubsteps 1000`, `-Pmin 0.0101` -- the
TIMs Workbench configuration.  This is a MECHANISM DEMONSTRATION, not a
bearing-capacity measurement: the mesh is far too coarse for a capacity and the
`-implex` reading hazard (ADR 92 section 8) applies to every number it prints.

WHAT IT SHOWS

  OLD binary (before WP-99):   every step "converges" and commits while the
                               companion counter of `implexRefusals` climbs by
                               (points x steps).  The load-settlement curve
                               keeps rising.  Nothing stops.

  NEW binary (WP-99):          the first capped commit refuses, nothing is
                               committed, the material latches, and the NEXT
                               step's update is refused -- `analyze()` returns
                               nonzero and the run stops.

Usage (from the worktree root):

    PYTHONPATH=<worktree>/dist/bin python3.12 \
        Ladruno_files/testbed/hypo_bearing/adr92_f7/strip_quad_implex.py out.csv

Optional second argument: number of settlement steps (default 40).
"""
import sys

import opensees as ops   # PYTHONPATH must point at a dist/bin

# ---------------------------------------------------------------------------
# CP1 / Gorini parameter set -- tests/test_ladruno_sanisand.py::_PARAMS
# ---------------------------------------------------------------------------
_PARAMS = [
    264.32, 0.3129, 0.6944, 1.33090, 0.71, 0.027, 0.83, 0.45, 101.0,
    0.005, 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0, 2.0,
]
_RHO = _PARAMS[17]          # 2.0 Mg/m3 with kPa/m/s units
_G = 9.81

# Geometry (half strip, symmetry at x = 0).  Coarse ON PURPOSE.
_B_HALF = 1.0               # half footing width [m]
_XMAX = 6.0
_ZMIN = -6.0
_H = 0.5                    # element size [m]

# The TIMs Workbench configuration, verbatim.
_MAXSUBSTEPS = 1000
_PMIN = 0.0101

_DS = 0.002                 # prescribed settlement increment [m/step]
_TOL = 1.0e-3               # relative force tolerance
_MAXITER = 40
_N_GRAV = 5


def _grid():
    nx = int(round(_XMAX / _H))
    nz = int(round(-_ZMIN / _H))
    nodes = {}
    tag = 1
    for i in range(nx + 1):
        for k in range(nz + 1):
            nodes[(i, k)] = tag
            ops.node(tag, i * _H, _ZMIN + k * _H)
            tag += 1
    return nx, nz, nodes


def build(mat_tag=9902, implex_control=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    nx, nz, nodes = _grid()

    opts = ["-Pmin", _PMIN, "-implex", "-maxSubsteps", _MAXSUBSTEPS]
    if implex_control:
        opts += ["-implexControl", 0.02]
    ops.nDMaterial("LadrunoSANISAND", mat_tag, *_PARAMS, *opts)

    eid = 1
    for i in range(nx):
        for k in range(nz):
            ops.element("LadrunoQuad", eid,
                        nodes[(i, k)], nodes[(i + 1, k)],
                        nodes[(i + 1, k + 1)], nodes[(i, k + 1)],
                        mat_tag,
                        "-thick", 1.0, "-type", "PlaneStrain",
                        "-formulation", "bbar",
                        "-rho", _RHO, "-body", 0.0, -_G)
            eid += 1

    # rollers on both vertical faces, pinned base -- ONE ops.fix per node
    # (a second fix on an already-constrained DOF is refused by the Domain).
    mask = {}
    for k in range(nz + 1):
        mask[nodes[(0, k)]] = [1, 0]
        mask[nodes[(nx, k)]] = [1, 0]
    for i in range(nx + 1):
        m = mask.setdefault(nodes[(i, 0)], [0, 0])
        m[1] = 1
    for n, m in mask.items():
        ops.fix(n, m[0], m[1])

    footing = [nodes[(i, nz)] for i in range(nx + 1) if i * _H <= _B_HALF + 1e-9]
    return nodes, footing, eid - 1


def _analysis(dlambda):
    # The tangent is UNSYMMETRIC (non-associated flow), and this mesh is tiny,
    # so a dense band solver is both adequate and free of setup surprises.
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-6, _MAXITER, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", dlambda)
    ops.analysis("Static")


def _refusals():
    r = list(ops.eleResponse(1, "material", 1, "implexRefusals"))
    while len(r) < 5:
        r.append(0.0)
    return r


def main(out_path, nsteps=40, implex_control=False):
    mat_tag = 9902
    nodes, footing, nele = build(mat_tag, implex_control)

    # --- gravity on the ELASTIC stage ------------------------------------
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)      # the -body self-weight is the load
    _analysis(1.0 / _N_GRAV)
    ops.updateMaterialStage("-material", mat_tag, "-stage", 0)
    for s in range(_N_GRAV):
        rc = ops.analyze(1)
        assert rc == 0, "gravity step %d failed (rc=%d)" % (s + 1, rc)
    ops.loadConst("-time", 0.0)
    ops.updateMaterialStage("-material", mat_tag, "-stage", 1)

    # --- prescribed settlement on the footing ----------------------------
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for n in footing:
        ops.sp(n, 2, -1.0)
    ops.integrator("LoadControl", _DS)

    rows = []
    stopped_at = None
    for s in range(nsteps):
        rc = ops.analyze(1)
        ops.reactions()
        q = sum(ops.nodeReaction(n, 2) for n in footing)
        w = -ops.nodeDisp(footing[0], 2)
        ref = _refusals()
        rows.append((s + 1, rc, w, -q, ref[0], ref[3], ref[4]))
        if rc != 0:
            stopped_at = s + 1
            break

    header = "step,rc,settlement_m,footing_load_kN_per_m,ref_total,ref_companion,commit_latched"
    with open(out_path, "w", newline="\n") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write("%d,%d,%.9g,%.9g,%d,%d,%d\n"
                     % (r[0], r[1], r[2], r[3], int(r[4]), int(r[5]), int(r[6])))

    print("build:", ops.ladrunoBuild().strip().splitlines()[0])
    print("elements: %d   footing nodes: %d" % (nele, len(footing)))
    print("wrote %s (%d rows)" % (out_path, len(rows)))
    if stopped_at is None:
        print("RESULT: ran all %d steps, every step committed. "
              "companion refusals = %d, latched = %d"
              % (nsteps, int(rows[-1][5]), int(rows[-1][6])))
    else:
        print("RESULT: STOPPED at step %d (rc != 0). "
              "companion refusals = %d, latched = %d"
              % (stopped_at, int(rows[-1][5]), int(rows[-1][6])))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "strip_quad_implex.csv"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    ctl = "--control" in sys.argv
    main(out, n, ctl)
