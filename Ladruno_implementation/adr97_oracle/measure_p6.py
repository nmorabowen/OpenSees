"""ADR-97 P6 -- Cerro-Lindo-scale measurement: Backward_Euler vs Closest_Point,
Secant/Continuum/Algorithmic tangents, on a real 3D bearing-capacity model.

WHAT THIS MEASURES
-------------------
A box of LadrunoBrick hexes (Mohr-Coulomb or MohrCoulombTensionCutoff,
non-associated psi < phi) is gravity-initialised (self-weight ramped, then
frozen with ``loadConst``) and then pushed toward bearing failure by a strip
footing pressure on the top-centre patch, ramped in LoadControl.  Five
configurations are compared:

    BE_Secant             Backward_Euler + Secant           (shipped default)
    BE_Continuum          Backward_Euler + Continuum
    CP_Continuum          Closest_Point  + Continuum
    CP_Algorithmic        Closest_Point  + Algorithmic       algorithm Newton
    CP_Algorithmic_Krylov Closest_Point  + Algorithmic       algorithm KrylovNewton

at two mesh scales: ``small`` (~2-5k DOF) and ``cerro`` (~22.6k DOF, the
ADR-80 hex8-class rung size -- see ``Ladruno_implementation/
80_ladruno_sp_imposition_strengthening_adr.md`` line 229, "All hex8-class
(~22.6 k DOF)").  A ``smoke`` scale (~150-300 DOF, 3 push steps) exists only
for ``tests/test_adr97_p6_measure_smoke.py``.

This is the P6 measurement input to ADR-97's D1 default-flip decision -- see
``Ladruno_implementation/97_ladruno_asdp_closest_point_adr.md`` D1 and the
Phases table's P6 row.

USAGE
-----
Single configuration, own process, JSON to stdout or --out:

    PYTHONPATH=<worktree>/dist/bin python3.12 measure_p6.py \\
        --family MC --scale small --config CP_Algorithmic --out run.json

Full sweep (spawns one FRESH child process per configuration -- required:
ASDPlasticMaterial3D's per-tag option maps are process-global statics, ADR-94
lesson, ADR-97 P1/P2 -- and aggregates into one JSON + a markdown table):

    PYTHONPATH=<worktree>/dist/bin python3.12 measure_p6.py --orchestrate \\
        --families MC,MCTC --scales small,cerro --out-dir <scratchpad>/p6

TRAPS OBSERVED (see ADR-97 docs / CLAUDE.md)
---------------------------------------------
* ``system FullGeneral`` crashes fully-prescribed rigs; not used here, but
  ``ProfileSPD`` is simply WRONG for these unsymmetric non-associated
  tangents -- always ``system UmfPack``.
* A ``Path`` time series returns 0 at its LAST point; gravity/footing ramps
  here use ``Linear`` series with ``loadConst`` freezing gravity, never Path.
* ``stdBrick`` swallows material return codes by design -- this driver only
  ever builds ``LadrunoBrick``.
* The material construction prints a banner/help page to stderr; harmless,
  discarded by the caller (child stdout/stderr is captured wholesale into the
  JSON's ``stderr_tail`` for debugging, never parsed for control flow).
* Per-tag static option maps: every (family, scale, config) triple below runs
  in ITS OWN ``python3.12`` child process, never reusing a material tag or an
  option map across configurations in one process.
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, os.pardir, os.pardir))

# ---------------------------------------------------------------------------
# material decks -- MC uses generic soil constants deliberately DIFFERENT from
# ADR-97 P2's oedometric trap (nu=0.25, phi=30 sit exactly on the compression
# meridian, K0 = nu/(1-nu) = 1/3 = sin(phi); a self-weight (K0) stress path
# with those exact numbers never yields).  MCTC reuses the Cerro-Lindo-like
# EDZ constants from tests/test_asdplastic_mctc.py / ADR-97 P2 (E=2.0e6 kPa,
# nu=0.3, c=100 kPa, phi=20, psi=5, T=24.7 kPa) -- a real geotechnical deck,
# not an invented one.
# ---------------------------------------------------------------------------
MC_E, MC_NU = 30000.0, 0.30
MC_PHI, MC_PSI, MC_C = 32.0, 8.0, 15.0
MC_GAMMA = 20.0  # kN/m^3 unit weight

TC_E, TC_NU = 2.0e6, 0.30
TC_C, TC_PHI, TC_PSI = 100.0, 20.0, 5.0
TC_T = 24.7
TC_GAMMA = 20.0

IV_NULL = "BackStress(NullHardeningTensorFunction):"

CONFIGS = {
    "BE_Secant":             dict(method="Backward_Euler", tangent="Secant",     algo="Newton"),
    "BE_Continuum":          dict(method="Backward_Euler", tangent="Continuum",  algo="Newton"),
    "CP_Continuum":          dict(method="Closest_Point",  tangent="Continuum",  algo="Newton"),
    "CP_Algorithmic":        dict(method="Closest_Point",  tangent="Algorithmic", algo="Newton"),
    "CP_Algorithmic_Krylov": dict(method="Closest_Point",  tangent="Algorithmic", algo="KrylovNewton"),
}
CONFIG_ORDER = list(CONFIGS.keys())

# scale -> (nx, ny, nz elements per axis, grav_steps, push_steps)
SCALES = {
    "smoke": dict(nx=3, ny=3, nz=2, grav_steps=2, push_steps=3),
    "small": dict(nx=10, ny=10, nz=10, grav_steps=3, push_steps=10),
    "cerro": dict(nx=19, ny=19, nz=19, grav_steps=3, push_steps=12),
}

LX, LY, LZ = 10.0, 10.0, 6.0   # metres
FOOTING_FRAC = 1.0 / 3.0        # central patch, plan fraction per axis
FOOTING_FACTOR_BY_FAMILY = {
    # fraction of the rough Terzaghi qult targeted by the end of the push
    # ramp -- comfortably into yield without demanding failure-point
    # convergence out of a coarse ramp; tuned empirically per family against
    # the small/smoke meshes.  MC's high-friction/low-cohesion deck yields
    # early (0.12*qult already drives 16-29 of 1000 GPs plastic); MCTC's
    # high-cohesion EDZ deck needs a much bigger push to leave the elastic
    # K0 state at all (0.12*qult left every GP elastic, measured).
    "MC": 0.12,
    "MCTC": 0.60,
}


def _terzaghi_qult(c, phi_deg, gamma, B):
    phi = math.radians(phi_deg)
    if phi < 1e-6:
        Nq, Nc, Ng = 1.0, 5.14, 0.0
    else:
        Nq = math.exp(math.pi * math.tan(phi)) * math.tan(math.pi / 4 + phi / 2) ** 2
        Nc = (Nq - 1.0) / math.tan(phi)
        Ng = 2.0 * (Nq + 1.0) * math.tan(phi)
    return c * Nc + 0.5 * gamma * B * Ng  # surface footing, q=0 overburden term


def _mat_deck(family):
    if family == "MC":
        return dict(E=MC_E, nu=MC_NU, phi=MC_PHI, psi=MC_PSI, c=MC_C, gamma=MC_GAMMA, T=None)
    if family == "MCTC":
        return dict(E=TC_E, nu=TC_NU, phi=TC_PHI, psi=TC_PSI, c=TC_C, gamma=TC_GAMMA, T=TC_T)
    raise ValueError("unknown family %r" % family)


def _make_material(ops, tag, family, method, tangent, niter=150):
    deck = _mat_deck(family)
    common = [
        "LinearIsotropic3D_EL", IV_NULL,
        "Begin_Model_Parameters",
        "YoungsModulus", deck["E"], "PoissonsRatio", deck["nu"],
        "MC_phi", deck["phi"], "MC_c", deck["c"], "MC_psi", deck["psi"],
        "MC_ds", 0.0,
    ]
    if family == "MCTC":
        common += ["TC_min_stress", deck["T"]]
    common += ["MassDensity", 0.0, "End_Model_Parameters",
               "Begin_Internal_Variables",
               "BackStress", 0., 0., 0., 0., 0., 0.,
               "End_Internal_Variables",
               "Begin_Integration_Options",
               "integration_method", method, "tangent_type", tangent,
               "n_max_iterations", int(niter),
               "End_Integration_Options"]
    if family == "MC":
        yf_pf = ["MohrCoulomb_YF", "MohrCoulomb_PF"]
    else:
        yf_pf = ["MohrCoulombTensionCutoff_YF", "MohrCoulombTensionCutoff_PF"]
    ops.nDMaterial("ASDPlasticMaterial3D", tag, *yf_pf, *common)


def node_tag(i, j, k, nx, ny):
    return 1 + i + j * (nx + 1) + k * (nx + 1) * (ny + 1)


def ele_tag(i, j, k, nx, ny):
    return 1 + i + j * nx + k * nx * ny


def build_and_run(ops, family, scale, config_name, out=None):
    t_start_all = time.perf_counter()
    scfg = SCALES[scale]
    nx, ny, nz = scfg["nx"], scfg["ny"], scfg["nz"]
    grav_steps, push_steps = scfg["grav_steps"], scfg["push_steps"]
    ccfg = CONFIGS[config_name]
    deck = _mat_deck(family)

    dx, dy, dz = LX / nx, LY / ny, LZ / nz

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    n_nodes = 0
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(node_tag(i, j, k, nx, ny), i * dx, j * dy, k * dz)
                n_nodes += 1

    _make_material(ops, 1, family, ccfg["method"], ccfg["tangent"])

    bz = -deck["gamma"]  # body force per unit volume, downward
    ele_tags = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                n1 = node_tag(i, j, k, nx, ny)
                n2 = node_tag(i + 1, j, k, nx, ny)
                n3 = node_tag(i + 1, j + 1, k, nx, ny)
                n4 = node_tag(i, j + 1, k, nx, ny)
                n5 = node_tag(i, j, k + 1, nx, ny)
                n6 = node_tag(i + 1, j, k + 1, nx, ny)
                n7 = node_tag(i + 1, j + 1, k + 1, nx, ny)
                n8 = node_tag(i, j + 1, k + 1, nx, ny)
                et = ele_tag(i, j, k, nx, ny)
                ops.element("LadrunoBrick", et, n1, n2, n3, n4, n5, n6, n7, n8,
                            1, "-b", 0.0, 0.0, bz)
                ele_tags.append(et)
    n_elements = len(ele_tags)

    # boundary conditions: base fully fixed, vertical side faces on rollers
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                tag = node_tag(i, j, k, nx, ny)
                if k == 0:
                    ops.fix(tag, 1, 1, 1)
                else:
                    fx = 1 if (i == 0 or i == nx) else 0
                    fy = 1 if (j == 0 or j == ny) else 0
                    if fx or fy:
                        ops.fix(tag, fx, fy, 0)

    n_dof = 3 * n_nodes  # before restraint removal; reported as model size

    # gravity pattern (self-weight, ramped then frozen)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", *ele_tags, "-type", "-selfWeight", 0.0, 0.0, 1.0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-8, 250, 0)
    ops.algorithm(ccfg["algo"])
    ops.integrator("LoadControl", 1.0 / grav_steps)
    ops.analysis("Static")

    grav = {"iters": [], "rc": [], "wall_time": 0.0, "first_fail_step": None}
    t0 = time.perf_counter()
    for s in range(grav_steps):
        rc = ops.analyze(1)
        grav["rc"].append(rc)
        if rc != 0:
            grav["first_fail_step"] = s
            break
        grav["iters"].append(ops.testIter())
    grav["wall_time"] = time.perf_counter() - t0

    result = {
        "family": family, "scale": scale, "config": config_name,
        "method": ccfg["method"], "tangent": ccfg["tangent"], "algorithm": ccfg["algo"],
        "nx": nx, "ny": ny, "nz": nz,
        "n_nodes": n_nodes, "n_elements": n_elements, "n_dof": n_dof,
        "grav_steps": grav_steps, "push_steps": push_steps,
        "grav": grav,
        "push": None,
        "monitor_ele": None,
        "footing_pressure_kPa": None,
        "plastic_gp_count": None,
        "plastic_gp_sampled": n_elements,
        "total_iters": sum(grav["iters"]),
        "total_wall_time": grav["wall_time"],
        "error": None,
    }

    if grav["first_fail_step"] is not None:
        result["total_wall_time"] = time.perf_counter() - t_start_all
        return result

    ops.loadConst("-time", 0.0)

    # footing patch: central FOOTING_FRAC of plan, on the top surface k=nz
    i_lo = max(1, round(nx * (0.5 - FOOTING_FRAC / 2)))
    i_hi = min(nx - 1, round(nx * (0.5 + FOOTING_FRAC / 2)))
    j_lo = max(1, round(ny * (0.5 - FOOTING_FRAC / 2)))
    j_hi = min(ny - 1, round(ny * (0.5 + FOOTING_FRAC / 2)))
    if i_hi <= i_lo:
        i_hi = i_lo + 1
    if j_hi <= j_lo:
        j_hi = j_lo + 1

    footing_nodes = [node_tag(i, j, nz, nx, ny)
                     for i in range(i_lo, i_hi + 1) for j in range(j_lo, j_hi + 1)]
    Bx = (i_hi - i_lo) * dx
    By = (j_hi - j_lo) * dy
    qult = _terzaghi_qult(deck["c"], deck["phi"], deck["gamma"], min(Bx, By))
    q_target = FOOTING_FACTOR_BY_FAMILY[family] * qult
    result["footing_pressure_kPa"] = q_target
    trib_area = (Bx * By) / len(footing_nodes)

    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for tag in footing_nodes:
        ops.load(tag, 0.0, 0.0, -q_target * trib_area)

    ops.test("NormDispIncr", 1.0e-8, 250, 0)
    ops.algorithm(ccfg["algo"])
    ops.integrator("LoadControl", 1.0 / push_steps)

    monitor_i = (i_lo + i_hi) // 2
    monitor_j = (j_lo + j_hi) // 2
    monitor_ele = ele_tag(monitor_i, monitor_j, nz - 1, nx, ny)
    monitor_node = node_tag(monitor_i, monitor_j, nz, nx, ny)
    result["monitor_ele"] = monitor_ele
    result["monitor_node"] = monitor_node

    push = {"iters": [], "rc": [], "disp_z": [], "monitor_stress": [],
            "wall_time": 0.0, "first_fail_step": None}
    t0 = time.perf_counter()
    for s in range(push_steps):
        rc = ops.analyze(1)
        push["rc"].append(rc)
        if rc != 0:
            push["first_fail_step"] = s
            break
        push["iters"].append(ops.testIter())
        push["disp_z"].append(ops.nodeDisp(monitor_node, 3))
        st = list(ops.eleResponse(monitor_ele, "stresses"))
        # average of the 8 GPs (48 = 8 x 6)
        avg = [sum(st[g * 6 + c] for g in range(8)) / 8.0 for c in range(6)]
        push["monitor_stress"].append(avg)
    push["wall_time"] = time.perf_counter() - t0
    result["push"] = push
    result["total_iters"] = sum(grav["iters"]) + sum(push["iters"])
    result["total_wall_time"] = time.perf_counter() - t_start_all

    # plastic-zone size: GP-1 sample of eqpstrain across every element (cheap,
    # one call per element instead of 8; documented approximation)
    try:
        n_plastic = 0
        for et in ele_tags:
            v = ops.eleResponse(et, "material", 1, "eqpstrain")
            if v and abs(v[0]) > 1.0e-10:
                n_plastic += 1
        result["plastic_gp_count"] = n_plastic
    except Exception as exc:  # pragma: no cover - diagnostic only
        result["plastic_gp_count"] = None
        result["plastic_zone_error"] = repr(exc)

    return result


def _child_main():
    p = argparse.ArgumentParser()
    p.add_argument("--family", required=True, choices=["MC", "MCTC"])
    p.add_argument("--scale", required=True, choices=list(SCALES.keys()))
    p.add_argument("--config", required=True, choices=CONFIG_ORDER)
    p.add_argument("--push-steps", type=int, default=None)
    p.add_argument("--grav-steps", type=int, default=None)
    p.add_argument("--out", default=None)
    args = p.parse_args()

    if args.push_steps is not None or args.grav_steps is not None:
        s = dict(SCALES[args.scale])
        if args.push_steps is not None:
            s["push_steps"] = args.push_steps
        if args.grav_steps is not None:
            s["grav_steps"] = args.grav_steps
        SCALES[args.scale] = s

    os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
    import opensees as ops  # noqa: E402  (import after env is set)

    build_hash = None
    try:
        build_hash = ops.ladrunoBuild()
    except Exception:
        pass

    try:
        result = build_and_run(ops, args.family, args.scale, args.config)
        result["build"] = build_hash
    except Exception as exc:
        result = {"family": args.family, "scale": args.scale, "config": args.config,
                   "error": repr(exc), "build": build_hash}

    text = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text)
    else:
        print(text)


# ---------------------------------------------------------------------------
# orchestrator: spawn one fresh child process per (family, scale, config)
# ---------------------------------------------------------------------------
def _orchestrate():
    p = argparse.ArgumentParser()
    p.add_argument("--families", default="MC,MCTC")
    p.add_argument("--scales", default="small,cerro")
    p.add_argument("--configs", default=",".join(CONFIG_ORDER))
    p.add_argument("--out-dir", required=True)
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--dist-dir", default=os.path.join(REPO_ROOT, "dist", "bin"))
    p.add_argument("--timeout", type=int, default=900)
    args = p.parse_args()

    families = args.families.split(",")
    scales = args.scales.split(",")
    configs = args.configs.split(",")
    os.makedirs(args.out_dir, exist_ok=True)

    env = dict(os.environ)
    env["PYTHONPATH"] = args.dist_dir + os.pathsep + env.get("PYTHONPATH", "")
    env["PATH"] = args.dist_dir + os.pathsep + env.get("PATH", "")
    env["LADRUNO_OPENSEES_QUIET"] = "1"

    all_results = []
    for family in families:
        for scale in scales:
            for config in configs:
                out_path = os.path.join(args.out_dir, "%s_%s_%s.json" % (family, scale, config))
                cmd = [args.python, os.path.abspath(__file__),
                       "--family", family, "--scale", scale, "--config", config,
                       "--out", out_path]
                sys.stderr.write("RUN %s\n" % " ".join(cmd))
                sys.stderr.flush()
                t0 = time.perf_counter()
                proc = subprocess.run(cmd, capture_output=True, text=True,
                                       env=env, timeout=args.timeout,
                                       stdin=subprocess.DEVNULL)
                wall = time.perf_counter() - t0
                if proc.returncode != 0 or not os.path.exists(out_path):
                    result = {"family": family, "scale": scale, "config": config,
                              "error": "child exit %d" % proc.returncode,
                              "stderr_tail": proc.stderr[-4000:],
                              "stdout_tail": proc.stdout[-2000:]}
                else:
                    with open(out_path) as fh:
                        result = json.load(fh)
                result["child_wall_time"] = wall
                all_results.append(result)
                sys.stderr.write("  -> done in %.1fs (child rc=%d)\n" % (wall, proc.returncode))
                sys.stderr.flush()

    combined = os.path.join(args.out_dir, "measure_p6_all.json")
    with open(combined, "w") as fh:
        json.dump(all_results, fh, indent=2)
    print("wrote", combined)
    return all_results


if __name__ == "__main__":
    if "--orchestrate" in sys.argv:
        sys.argv.remove("--orchestrate")
        _orchestrate()
    else:
        _child_main()
