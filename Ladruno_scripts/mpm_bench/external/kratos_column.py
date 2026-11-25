#!/usr/bin/env python
"""Matched elastic explicit 3D column in Kratos Multiphysics (MPMApplication, BSD), WP-253.

Case = Ladruno_scripts/mpm_m0b/models.py `column_commands(N, "elastic")`: quarter column
a x a x 2a (a = 1 m), cell h = a/m, particles per cell 2 per axis (8), E = 1e6 Pa, nu = 0.3,
rho = 2000 kg/m3, g = 9.81 down, dt = 0.25 dt_crit, symmetry planes x = 0 and y = 0,
floor at z = 0. N = 8 * 9*9*18 = 11,664 (m = 9) or 93,312 (m = 18).

Usage:  python kratos_column.py --N 10000 --threads 1 [--steps 20] [--warmup 5] [--scheme musl]
Prints one JSON line (Kratos's teardown warnings follow it on stdout: grep '^{', or use
--json-out FILE). The timing window is the wall of the full per-step cycle
(AdvanceInTime + InitializeSolutionStep + Predict (search + P2G) + SolveSolutionStep
+ FinalizeSolutionStep); no output process, no mesh/particle generation (that is in Initialize,
before the clock). See README.md for what is and is not matched.
"""
import argparse, json, math, os, sys, tempfile, time

M_FOR_N = {1000: 4, 2000: 5, 10000: 9, 100000: 18}
A, RHO, G, E, NU, DTF = 1.0, 2000.0, 9.81, 1.0e6, 0.3, 0.25


def write_mdpa(path, nodes, elems, ename, prop_id=1, submodel=None):
    with open(path, "w") as f:
        f.write("Begin Properties %d\nEnd Properties\n\nBegin Nodes\n" % prop_id)
        for i, (x, y, z) in enumerate(nodes, 1):
            f.write("%d %.12g %.12g %.12g\n" % (i, x, y, z))
        f.write("End Nodes\n\nBegin Elements %s\n" % ename)
        for i, c in enumerate(elems, 1):
            f.write("%d %d %s\n" % (i, prop_id, " ".join(str(n) for n in c)))
        f.write("End Elements\n")
        if submodel:
            f.write("\nBegin SubModelPart %s\n Begin SubModelPartNodes\n" % submodel)
            for i in range(1, len(nodes) + 1):
                f.write(" %d\n" % i)
            f.write(" End SubModelPartNodes\n Begin SubModelPartElements\n")
            for i in range(1, len(elems) + 1):
                f.write(" %d\n" % i)
            f.write(" End SubModelPartElements\nEnd SubModelPart\n")


def hex_mesh(nx, ny, nz, h):
    nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)
    nodes = [(i * h, j * h, k * h) for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)]
    el = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                el.append((nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                           nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)))
    return nodes, el


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--steps", type=int, default=20)
    ap.add_argument("--warmup", type=int, default=5)
    ap.add_argument("--scheme", default="musl", choices=["musl", "usl", "usf", "central_difference"])
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--json-out", default=None, help="also write the result JSON to this file")
    ap.add_argument("--law", default="HyperElasticNeoHookean3DLaw")
    a = ap.parse_args()
    # thread env must be set before Kratos loads
    os.environ["OMP_NUM_THREADS"] = str(a.threads)
    os.environ["MKL_NUM_THREADS"] = str(a.threads)
    m = M_FOR_N[min(M_FOR_N, key=lambda k: abs(k - a.N))]
    h = A / m
    wd = a.workdir or tempfile.mkdtemp(prefix="kratos_col_")
    os.makedirs(wd, exist_ok=True)
    os.chdir(wd)

    import KratosMultiphysics as KM
    import KratosMultiphysics.MPMApplication as KMPM
    from KratosMultiphysics.MPMApplication.mpm_analysis import MpmAnalysis
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    try:
        KM.ParallelUtilities.SetNumThreads(a.threads)
    except Exception:
        pass

    gnodes, gel = hex_mesh(4 * m, 4 * m, 2 * m + 2, h)      # fork grid: 4m x 4m x (2m+2) cells
    write_mdpa("Grid.mdpa", gnodes, gel, "Element3D8N")
    bnodes, bel = hex_mesh(m, m, 2 * m, h)                  # body: m x m x 2m cells, 8 MP each
    write_mdpa("Body.mdpa", bnodes, bel, "Element3D8N", submodel="Parts_Material_domain")
    n_particles = 8 * len(bel)

    # dt_crit = h / c, c = sqrt(M / rho), M = constrained modulus; dt = 0.25 dt_crit
    M = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
    c = math.sqrt(M / RHO)
    dt_c = h / c
    dt = DTF * dt_c

    mats = {"properties": [{"model_part_name": "Initial_MPM_Material.Parts_Material_domain", "properties_id": 1,
            "Material": {"constitutive_law": {"name": a.law},
                         "Variables": {"DENSITY": RHO, "YOUNG_MODULUS": E, "POISSON_RATIO": NU,
                                       "MATERIAL_POINTS_PER_ELEMENT": 8}, "Tables": {}}}]}
    json.dump(mats, open("Materials.json", "w"))
    stress = {"musl": "musl", "usl": "usl", "usf": "usf", "central_difference": "usf"}[a.scheme]
    solver = {
        "solver_type": "Dynamic", "model_part_name": "MPM_Material", "domain_size": 3, "echo_level": 0,
        "analysis_type": "linear", "time_integration_method": "explicit",
        "scheme_type": "central_difference" if a.scheme == "central_difference" else "forward_euler",
        "stress_update": stress,
        "time_stepping": {"time_step": dt},
        "model_import_settings": {"input_type": "mdpa", "input_filename": "Body"},
        "material_import_settings": {"materials_filename": "Materials.json"},
        "grid_model_import_settings": {"input_type": "mdpa", "input_filename": "Grid"},
        "pressure_dofs": False, "compressible": True,
        "problem_domain_sub_model_part_list": ["Parts_Material_domain"],
        "processes_sub_model_part_list": [], "auxiliary_variables_list": [],
        "linear_solver_settings": {"solver_type": "skyline_lu_factorization"},
    }
    params = KM.Parameters(json.dumps({
        "problem_data": {"problem_name": "col", "parallel_type": "OpenMP", "echo_level": 0,
                         "start_time": 0.0, "end_time": 1.0e9},
        "solver_settings": solver, "processes": {}, "output_processes": {}}))
    model = KM.Model()
    an = MpmAnalysis(model, params)
    t0 = time.perf_counter()
    an.Initialize()
    t_init = time.perf_counter() - t0
    grid = model["Background_Grid"]
    mp = model["MPM_Material"]
    # symmetry planes and floor (grid-node Dirichlet; see README: floor is frictionless)
    eps = 1e-9 * h
    for n in grid.Nodes:
        if n.X < eps:
            n.Fix(KM.DISPLACEMENT_X); n.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
        if n.Y < eps:
            n.Fix(KM.DISPLACEMENT_Y); n.SetSolutionStepValue(KM.DISPLACEMENT_Y, 0.0)
        if n.Z < eps:
            n.Fix(KM.DISPLACEMENT_Z); n.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0.0)
    gvec = KM.Vector([0.0, 0.0, -G])
    for el in mp.Elements:
        el.SetValuesOnIntegrationPoints(KMPM.MP_VOLUME_ACCELERATION, [gvec], mp.ProcessInfo)
        el.SetValuesOnIntegrationPoints(KMPM.MP_ACCELERATION, [gvec], mp.ProcessInfo)
    n_mp0 = mp.NumberOfElements()

    solver_obj = an._GetSolver()

    def one_step():
        an.time = solver_obj.AdvanceInTime(an.time)
        an.InitializeSolutionStep()
        solver_obj.Predict()
        solver_obj.SolveSolutionStep()
        an.FinalizeSolutionStep()

    an.time = 0.0
    for _ in range(a.warmup):
        one_step()
    walls = []
    for _ in range(a.steps):
        t = time.perf_counter(); one_step(); walls.append(time.perf_counter() - t)
    W = sum(walls)
    n_mp1 = mp.NumberOfElements()
    ws = sorted(walls)
    # sanity (after the clock): mean MP vertical velocity and lowest/highest MP z
    zs, vz = [], []
    for el in mp.Elements:
        zs.append(el.CalculateOnIntegrationPoints(KMPM.MP_COORD, mp.ProcessInfo)[0][2])
        vz.append(el.CalculateOnIntegrationPoints(KMPM.MP_VELOCITY, mp.ProcessInfo)[0][2])
    out = {"mp_z_min": min(zs), "mp_z_max": max(zs), "mp_vz_mean": sum(vz) / len(vz),
           "max_abs_vz": max(abs(v) for v in vz), "law": a.law,
           "code": "Kratos MPMApplication", "kratos_version": getattr(KM, "KRATOS_VERSION", "") or "10.4.4",
           "python": sys.version.split()[0], "scheme": a.scheme, "N_requested": a.N, "N": n_particles,
           "particles_after": n_mp1, "particles_start": n_mp0, "m_cells": m, "h": h, "dt": dt, "dt_crit": dt_c,
           "threads": a.threads, "kratos_threads": KM.ParallelUtilities.GetNumThreads(),
           "warmup": a.warmup, "steps": a.steps, "init_s": t_init, "window_s": W,
           "step_median_s": ws[len(ws) // 2],
           "throughput_particle_steps_per_s": n_mp0 * a.steps / W,
           "throughput_from_median": n_mp0 / ws[len(ws) // 2]}
    # Kratos writes its teardown warnings to stdout AFTER this line: read the line that starts
    # with "{" (or --json-out), never the last line
    print(json.dumps(out), flush=True)
    if a.json_out:
        with open(a.json_out, "w") as f:
            f.write(json.dumps(out) + "\n")


if __name__ == "__main__":
    main()
