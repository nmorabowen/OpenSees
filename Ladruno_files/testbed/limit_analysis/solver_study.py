"""WP-237d: P0 solver and certification study driver.

Every job is one (case, mesh, formulation, solver) run in its own worker process under a supervisor that
enforces a wall-clock limit and a resident-memory guard and records the peak working set. The worker builds
the program with the frozen reference (``tests/limit_reference.py``), solves it, writes problem, result and
interior-field files in the frozen formats, then runs the independent checker (``tests/limit_certificate.py``)
on the files it wrote and keeps each certificate with its ``.certificate.method.json``. Gate evidence is read
only through ``limit_certificate.read_gate_certificate``.

Nothing here edits the reference, the checker or the contract. Every number a run produces is an estimate
until the checker certifies it.

Commands::

    python solver_study.py job    --case C --mesh M --formulation F --solver S --data DIR [--time-limit s] [--mem-gb g]
    python solver_study.py queue  --file JOBS.txt --data DIR           one job per line: case mesh formulation solver
    python solver_study.py table  --data DIR --out SUMMARY.json [--md SUMMARY.md]
    python solver_study.py sizes  --out FILE                           program sizes without solving (assemble only)

Solvers: ``clarabel-qdldl``, ``clarabel-faer``, ``scs``, ``qoco`` (reported ``not_installed`` when absent).
Large data (problems, results, certificates) goes in ``--data``, outside git; the small row files are the
record.
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
TESTS = HERE.parents[2] / "tests"
for path in (str(TESTS),):
    if path not in sys.path:
        sys.path.insert(0, path)

ROW_FORMAT = "ladruno-limit-solver-study-row"
ROW_VERSION = 0

# Named budgets. Not certificate tolerances. Calibration: the brief caps a run at 30 minutes and 20 GB; this
# workstation has about 3 GB of free memory when the study starts, so the guard is far below the cap.
TIME_LIMIT_S = 1500.0          # solver wall-clock stop inside the worker; the supervisor kills at +300 s
SUPERVISOR_SLACK_S = 300.0
MEMORY_GUARD_GB = 2.6          # resident-set guard of the supervisor
POLL_S = 0.2

PRESETS = {"nodyn": {"dynamic_regularization_enable": False}}   # Clarabel option variants named in queue files
SOLVERS = ("clarabel-qdldl", "clarabel-faer", "scs", "qoco")
SCS_MAX_ITERS = 100_000
SCS_EPS = 1e-6                 # SCS default eps_abs, eps_rel at 1e-4 in 3.2.x; a cross-check, never a certificate


def stem_of(case, mesh, formulation, solver, tag=""):
    return f"{case}.m{mesh}.q.{formulation}.{solver}" + (f".{tag}" if tag else "")


# ---------------------------------------------------------------------------------------------------
# Solvers: one canonical program  min c^T x,  A x + s = b,  s in K  (interface 15.4)
# ---------------------------------------------------------------------------------------------------

def peak_gb():
    import psutil
    info = psutil.Process().memory_info()
    return getattr(info, "peak_wset", info.rss) / 2 ** 30


def solve_clarabel(program, method, time_limit, overrides=None):
    import clarabel
    import numpy as np
    import scipy.sparse as sp
    import limit_reference as ref
    cones = ref._clarabel_cones(clarabel, program.n_zero, program.cones)
    n = program.A.shape[1]
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.time_limit = float(time_limit)
    settings.direct_solve_method = method
    for key, value in (overrides or {}).items():
        setattr(settings, key, value)
    solver = clarabel.DefaultSolver(sp.csc_matrix((n, n)), np.asarray(program.c, float), sp.csc_matrix(program.A),
                                    np.asarray(program.b, float), cones, settings)
    start = time.perf_counter()
    solution = solver.solve()
    wall = time.perf_counter() - start
    raw = str(solution.status).split(".")[-1]
    status = ref.CLARABEL_STATUS.get(raw, "numerical_error")
    info = {
        "status_raw": raw, "iterations": int(solution.iterations), "solve_time_s": float(solution.solve_time),
        "r_prim": float(solution.r_prim), "r_dual": float(solution.r_dual),
        "obj_primal": float(solution.obj_val), "obj_dual": float(solution.obj_val_dual),
        "rel_gap": abs(float(solution.obj_val) - float(solution.obj_val_dual))
        / max(1.0, abs(float(solution.obj_val)), abs(float(solution.obj_val_dual))),
    }
    record = {"name": "clarabel", "version": str(clarabel.__version__),
              "settings": {"defaults": True, "direct_solve_method": method, "overrides": dict(overrides or {}), "time_limit": float(time_limit),
                           "verbose": False, "tol_gap_abs": float(settings.tol_gap_abs),
                           "tol_gap_rel": float(settings.tol_gap_rel), "tol_feas": float(settings.tol_feas),
                           "max_iter": int(settings.max_iter), "status_raw": raw}}
    return status, np.asarray(solution.x, dtype=float), wall, record, info


SCS_STATUS = {"solved": "solved", "solved_inaccurate": "almost_solved", "infeasible": "primal_infeasible",
              "unbounded": "dual_infeasible", "infeasible_inaccurate": "numerical_error",
              "unbounded_inaccurate": "numerical_error", "failure": "numerical_error"}


def solve_scs(program, time_limit):
    import numpy as np
    import scipy.sparse as sp
    import scs
    q = []
    for cone in program.cones:
        q.extend([int(cone["dim"])] * int(cone["count"]))
    data = {"A": sp.csc_matrix(program.A), "b": np.asarray(program.b, float), "c": np.asarray(program.c, float)}
    cone = {"z": int(program.n_zero), "q": q}
    start = time.perf_counter()
    solver = scs.SCS(data, cone, max_iters=SCS_MAX_ITERS, eps_abs=SCS_EPS, eps_rel=SCS_EPS,
                     time_limit_secs=float(time_limit), verbose=False)
    solution = solver.solve()
    wall = time.perf_counter() - start
    info = solution["info"]
    raw = str(info["status"])
    status = SCS_STATUS.get(raw, "numerical_error")
    if raw == "solved (inaccurate - reached max_iters)" or "max_iters" in raw:
        status = "iteration_limit"
    if "time" in raw:
        status = "time_limit"
    x = np.asarray(solution["x"], dtype=float)
    out = {"status_raw": raw, "iterations": int(info["iter"]), "solve_time_s": float(info["solve_time"]) / 1000.0,
           "r_prim": float(info["res_pri"]), "r_dual": float(info["res_dual"]),
           "obj_primal": float(info["pobj"]), "obj_dual": float(info["dobj"]), "rel_gap": float(info["gap"])}
    record = {"name": "scs", "version": str(scs.__version__),
              "settings": {"eps_abs": SCS_EPS, "eps_rel": SCS_EPS, "max_iters": SCS_MAX_ITERS,
                           "time_limit_secs": float(time_limit), "status_raw": raw}}
    return status, x, wall, record, out


def solve_qoco(program, time_limit):
    """QOCO takes ``A_eq x = b_eq``, ``G x <=_K h`` with nonnegative and second-order cones only. The frozen
    programs have a zero cone and second-order cones, so they fit; no PSD cone exists in version 0."""
    try:
        import qoco
    except ImportError:
        return "not_run", None, 0.0, None, {"status_raw": "qoco_not_installed"}
    import numpy as np
    import scipy.sparse as sp
    nz = int(program.n_zero)
    A = sp.csc_matrix(program.A)
    Aeq, G = A[:nz], A[nz:]
    q = []
    for cone in program.cones:
        q.extend([int(cone["dim"])] * int(cone["count"]))
    model = qoco.QOCO()
    n = A.shape[1]
    start = time.perf_counter()
    model.setup(n, G.shape[0], Aeq.shape[0], sp.csc_matrix((n, n)), np.asarray(program.c, float), sp.csc_matrix(Aeq),
                np.asarray(program.b[:nz], float), sp.csc_matrix(G), np.asarray(program.b[nz:], float), 0, len(q),
                np.asarray(q, dtype=np.int32), {"verbose": 0})
    solution = model.solve()
    wall = time.perf_counter() - start
    raw = str(solution.status)
    status = {"QOCO_SOLVED": "solved", "QOCO_SOLVED_INACCURATE": "almost_solved",
              "QOCO_MAX_ITER": "iteration_limit"}.get(raw, "numerical_error")
    out = {"status_raw": raw, "iterations": int(solution.iters), "solve_time_s": float(solution.solve_time_sec),
           "obj_primal": float(solution.obj)}
    record = {"name": "qoco", "version": str(getattr(qoco, "__version__", "?")), "settings": {"verbose": 0}}
    return status, np.asarray(solution.x, dtype=float), wall, record, out


def solve_program(program, solver, time_limit, overrides=None):
    if solver == "clarabel-qdldl":
        return solve_clarabel(program, "qdldl", time_limit, overrides)
    if solver == "clarabel-faer":
        return solve_clarabel(program, "faer", time_limit, overrides)
    if solver == "scs":
        return solve_scs(program, time_limit)
    if solver == "qoco":
        return solve_qoco(program, time_limit)
    raise SystemExit(f"unknown solver {solver}")


# ---------------------------------------------------------------------------------------------------
# Worker: one job
# ---------------------------------------------------------------------------------------------------

def patch_interior_solves(overrides):
    """Give the reference's interior-field solves the same Clarabel options as the main solve.

    The reference calls ``_clarabel_solve`` with library defaults. This replaces that one function in the worker
    process only; no file of the reference changes. It is what the study's ``nodyn`` variant needs: the interior
    program of a frictional Upper is itself a conic program that stalls under dynamic regularisation.
    """
    import clarabel
    import numpy as np
    import scipy.sparse as sp
    import limit_reference as ref

    def patched(matrix, b, c, cones, time_limit):
        n = matrix.shape[1]
        settings = clarabel.DefaultSettings()
        settings.verbose = False
        settings.time_limit = float(time_limit)
        for key, value in overrides.items():
            setattr(settings, key, value)
        start = time.perf_counter()
        solution = clarabel.DefaultSolver(sp.csc_matrix((n, n)), np.asarray(c, float), sp.csc_matrix(matrix),
                                          np.asarray(b, float), cones, settings).solve()
        raw = str(solution.status).split(".")[-1]
        return (ref.CLARABEL_STATUS.get(raw, "numerical_error"), np.asarray(solution.x, dtype=float),
                time.perf_counter() - start, {"status_raw": raw})
    ref._clarabel_solve = patched


def needs_lower_interior(geo, program):
    """The checker takes a zero interior field for a homogeneous problem with every ``k > 0`` (interface 17.3)."""
    import numpy as np
    import limit_reference as ref
    f0 = ref.query_tables(program.problem, program.query)[1]
    homogeneous = not np.any(geo.b0) and not np.any(np.asarray(geo.arrays["t0"])) and not np.any(f0)
    return not (homogeneous and np.all(geo.k > 0.0))


def worker(case, mesh, formulation, solver, data, time_limit, overrides=None, tag="", solve_only=False):
    import numpy as np
    import limit_certificate as lc
    import limit_reference as ref
    data = Path(data)
    data.mkdir(parents=True, exist_ok=True)
    row = {"format": ROW_FORMAT, "format_version": ROW_VERSION, "case": case, "mesh": mesh,
           "formulation": formulation, "solver": solver, "phases": {}, "peak_gb": {}}
    if overrides:
        patch_interior_solves(overrides)
    stem = stem_of(case, mesh, formulation, solver, tag)
    row["tag"], row["overrides"] = tag, dict(overrides or {})

    def phase(name, start):
        row["phases"][name] = round(time.perf_counter() - start, 3)
        row["peak_gb"][name] = round(peak_gb(), 3)

    def save():
        (data / f"{stem}.row.json").write_text(json.dumps(row, indent=1, sort_keys=True) + "\n", encoding="utf-8")

    try:
        t = time.perf_counter()
        problem, arrays = ref.case_problem(case, mesh)
        ref.validate(problem, arrays)
        geo = ref.Geometry(problem, arrays)
        row["n_tets"] = int(geo.E)
        phase("build", t)
        t = time.perf_counter()
        program = ref.assemble(problem, arrays, formulation, "q", geometry=geo, validated=True)
        row.update(n_rows=int(program.A.shape[0]), n_vars=int(program.A.shape[1]), nnz=int(program.A.nnz),
                   n_zero=int(program.n_zero), bound=program.bound)
        phase("assemble", t)
        save()
        ref.export_problem(data, problem, arrays)

        t = time.perf_counter()
        status, x, wall, record, info = solve_program(program, solver, time_limit, overrides)
        phase("solve", t)
        row.update(solver_status=status, solver_wall_s=round(wall, 3), solver_record=record, solver_info=info)
        if x is not None and status in ("solved", "almost_solved"):
            row["zero_residual_inf"] = float(np.max(np.abs(program.zero_residual(x)), initial=0.0))
        save()
        if solve_only and x is not None:
            row["value_estimate"] = float(program.value(x)) if status in ("solved", "almost_solved") else None
            return
        if x is None or status not in ("solved", "almost_solved"):
            row["certificate"] = {"classification": "no_value", "note": f"solver_status {status}"}
            save()
            return

        extra = None
        t = time.perf_counter()
        frictional = bool(np.any(geo.alpha > 0.0))
        try:
            if program.bound == "upper" and frictional:
                extra = ref.upper_interior_field(program, time_limit)[:2]
            elif program.bound == "lower" and needs_lower_interior(geo, program):
                found = ref.lower_interior_field(program, time_limit)
                extra = None if found is None else found[:2]
            row["interior"] = "supplied" if extra is not None else "none"
        except Exception as error:   # an interior-field failure is a result, not a crash
            row["interior"] = f"failed: {type(error).__name__}: {error}"
        phase("interior", t)
        save()

        doc, result_arrays = ref.export_result(data, program, status, x, record, wall, extra)
        row["value_estimate"] = doc["value"]
        row["objective"] = doc["objective"]
        row["result_status"] = doc["solver_status"]
        del program, x
        save()
        if doc["value"] is None:
            row["certificate"] = {"classification": "no_value", "note": f"result status {doc['solver_status']}"}
            return

        t = time.perf_counter()
        stem_ref = ref.record_stem(problem, "q", formulation)
        # The checker reads the files the reference wrote, as an independent consumer would.
        cproblem, carrays = lc.load_document(data / f"{ref.problem_stem(problem)}.problem.json")
        cresult, cfields = lc.load_document(data / f"{stem_ref}.result.json")
        try:
            outcome = lc.check(cproblem, carrays, cresult, cfields)
        except lc.CapacityRefusal as error:
            row["certificate"] = {"classification": "capacity_refusal", "note": str(error)[:400]}
            phase("check", t)
            return
        except lc.Refusal as error:
            row["certificate"] = {"classification": "refusal", "note": f"{error}"[:400]}
            phase("check", t)
            return
        phase("check", t)
        cert = outcome.certificate
        row["certificate"] = {
            "classification": cert["classification"], "reasons": list(cert["reasons"]),
            "value_certified": cert["value_certified"], "floor": cert["floor"],
            "domain_claim": cert["domain_claim"], "checks": cert["checks"],
        }
        row["kappa_eq"] = cert["checks"]["kappa_eq"]
        row["repair_mu"] = cert["checks"]["repair_mu"]
        # One file name per solver so that solvers on one problem do not overwrite each other.
        directory = data / "certs" / (solver + (f".{tag}" if tag else ""))
        directory.mkdir(parents=True, exist_ok=True)
        if outcome.certificate is not None and outcome.method is not None:
            path = lc.write_certificate(directory, cresult, outcome)
            row["certificate_file"] = str(Path("certs") / directory.name / path.name)
            try:
                gate_cert, gate_method = lc.read_gate_certificate(path)
                row["gate_admissible"] = True
                row["interior_field"] = gate_method["interior_field"]
                row["kappa_method"] = gate_method["kappa_method"]
            except lc.Refusal as error:
                row["gate_admissible"] = False
                row["gate_refusal"] = str(error)[:300]
        # The result files are large; keep only the certificate directory for the gate record.
        (data / f"{stem_ref}.result.json").replace(directory / f"{stem_ref}.result.json")
        (data / f"{stem_ref}.result.npz").replace(directory / f"{stem_ref}.result.npz")
        row["peak_gb"]["final"] = round(peak_gb(), 3)
    except MemoryError as error:
        row["error"] = {"type": "MemoryError", "message": str(error)[:300], "trace": traceback.format_exc()[-1500:]}
    except Exception as error:
        row["error"] = {"type": type(error).__name__, "message": str(error)[:300], "trace": traceback.format_exc()[-1500:]}
    finally:
        row["peak_gb"]["final"] = round(peak_gb(), 3)
        save()


# ---------------------------------------------------------------------------------------------------
# Supervisor
# ---------------------------------------------------------------------------------------------------

def supervise(case, mesh, formulation, solver, data, time_limit=TIME_LIMIT_S, mem_gb=MEMORY_GUARD_GB, overrides=None, tag="", solve_only=False):
    import psutil
    data = Path(data)
    data.mkdir(parents=True, exist_ok=True)
    stem = stem_of(case, mesh, formulation, solver, tag)
    row_path = data / f"{stem}.row.json"
    if row_path.exists():
        row_path.unlink()
    command = [sys.executable, str(Path(__file__).resolve()), "worker", "--case", case, "--mesh", str(mesh),
               "--formulation", formulation, "--solver", solver, "--data", str(data), "--time-limit", str(time_limit),
               "--overrides", json.dumps(overrides or {}), "--tag", tag] + (["--solve-only"] if solve_only else [])
    env = dict(os.environ, LADRUNO_OPENSEES_QUIET="1", APEGMSH_QUIET="1", PYTHONUNBUFFERED="1")
    start = time.perf_counter()
    process = subprocess.Popen(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    handle = psutil.Process(process.pid)
    peak, outcome = 0.0, "completed"
    while process.poll() is None:
        try:
            rss = sum(item.memory_info().rss for item in [handle] + handle.children(recursive=True)) / 2 ** 30
        except psutil.Error:
            break
        peak = max(peak, rss)
        if rss > mem_gb:
            outcome = "memory_guard"
            process.kill()
            break
        if time.perf_counter() - start > time_limit + SUPERVISOR_SLACK_S:
            outcome = "timeout"
            process.kill()
            break
        time.sleep(POLL_S)
    try:
        tail = process.communicate(timeout=30)[0][-600:]
    except Exception:
        tail = ""
    wall = time.perf_counter() - start
    row = json.loads(row_path.read_text(encoding="utf-8")) if row_path.exists() else {
        "format": ROW_FORMAT, "format_version": ROW_VERSION, "case": case, "mesh": mesh,
        "formulation": formulation, "solver": solver, "phases": {}, "peak_gb": {}}
    row.update(outcome=outcome, supervisor_wall_s=round(wall, 2), supervisor_peak_rss_gb=round(peak, 3),
               exit_code=process.returncode, memory_guard_gb=mem_gb)
    if outcome == "completed" and process.returncode not in (0, None):
        row["outcome"] = "crash"
        row["crash_tail"] = tail
    if outcome != "completed":
        row.setdefault("certificate", {"classification": "no_value", "note": outcome})
    row_path.write_text(json.dumps(row, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    cert = row.get("certificate", {})
    print(f"{case} m{mesh} {formulation} {solver} {tag}: {row['outcome']} status={row.get('solver_status')} "
          f"it={row.get('solver_info', {}).get('iterations')} wall={wall:.1f}s peak={peak:.2f}GB "
          f"cert={cert.get('classification')} kappa={row.get('kappa_eq')}", flush=True)
    return row


def queue(path, data, time_limit, mem_gb):
    for line in Path(path).read_text().splitlines():
        line = line.split("#")[0].split()
        if len(line) in (4, 5):
            tag = line[4] if len(line) == 5 else ""
            supervise(line[0], int(line[1]), line[2], line[3], data, time_limit, mem_gb, PRESETS.get(tag), tag)


# ---------------------------------------------------------------------------------------------------
# Table and sizes
# ---------------------------------------------------------------------------------------------------

def table(data, out, md=None):
    """Summary of every row, with the best certified bracket per (case, mesh).

    Certified values are read only through ``read_gate_certificate``. The best Lower is the largest certified Lower
    of any solver run, the best Upper the smallest certified Upper of one formulation; each is a certified bound on
    its own, so the pair is a bracket as the checker's ``bracket`` decides, never as arithmetic here.
    """
    import limit_certificate as lc
    data = Path(data)
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted(data.glob("*.row.json"))]
    gate = {}
    for index, row in enumerate(rows):
        row["label"] = row["solver"] + (f".{row['tag']}" if row.get("tag") else "")
        file = row.get("certificate_file")
        row["gate_value"] = None
        if file:
            try:
                certificate, method = lc.read_gate_certificate(data / file)
                row["gate_value"], row["gate_floor"] = certificate["value_certified"], certificate["floor"]
                row["kappa_eq"] = certificate["checks"]["kappa_eq"]
                gate[index] = (certificate, method)
            except lc.Refusal as error:
                row["gate_refusal"] = str(error)[:200]
    best = {}
    for index, (certificate, method) in gate.items():
        row = rows[index]
        key = (row["case"], row["mesh"])
        kind = "lower" if row["bound"] == "lower" else ("jump" if row["formulation"] == lc.JUMP_ID else "control")
        entry = best.setdefault(key, {})
        value = certificate["value_certified"]
        if kind not in entry or (kind == "lower" and value > entry[kind][0]) or (kind != "lower" and value < entry[kind][0]):
            entry[kind] = (value, certificate, method, row["label"])
    brackets = []
    for (case, mesh), entry in sorted(best.items()):
        for kind in ("jump", "control"):
            if "lower" in entry and kind in entry:
                low, up = entry["lower"], entry[kind]
                item = {"case": case, "mesh": mesh, "upper_formulation": up[1]["formulation"], "lower": low[0],
                        "upper": up[0], "lower_run": low[3], "upper_run": up[3]}
                try:
                    record = lc.bracket(low[1], up[1], low[2], up[2])
                    item.update(status=record["status"], gap=record["gap"], domain_claim=record["domain_claim"])
                except lc.Refusal as error:
                    item["status"] = f"refused: {error}"[:120]
                brackets.append(item)
        for kind in ("lower", "jump", "control"):
            if kind in entry and not any(b["case"] == case and b["mesh"] == mesh for b in brackets):
                brackets.append({"case": case, "mesh": mesh, "single": kind, "value": entry[kind][0], "run": entry[kind][3]})
    summary = {"format": "ladruno-limit-solver-study-summary", "format_version": 0, "rows": rows, "brackets": brackets}
    Path(out).write_text(json.dumps(summary, indent=1, sort_keys=True) + chr(10), encoding="utf-8")
    if md:
        lines = ["| case | m | formulation | solver | outcome | status | it | solve s | peak GB | cert | value | kappa_eq |",
                 "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in rows:
            info = r.get("solver_info", {})
            lines.append("| {case} | {mesh} | {f} | {solver} | {o} | {s} | {it} | {w} | {pk} | {c} | {v} | {k} |".format(
                case=r["case"], mesh=r["mesh"], f=r["formulation"].split(".")[1], solver=r["label"],
                o=r.get("outcome"), s=r.get("solver_status"), it=info.get("iterations"),
                w=r.get("solver_wall_s"), pk=r.get("supervisor_peak_rss_gb"),
                c=r.get("certificate", {}).get("classification"),
                v=None if r.get("gate_value") is None else f"{r['gate_value']:.6g}", k=r.get("kappa_eq")))
        Path(md).write_text(chr(10).join(lines) + chr(10), encoding="utf-8")
    return summary


def sizes(out, cases=None, meshes=None):
    import limit_reference as ref
    import test_limit_contract as contract
    names = cases or [n for n, c in contract.parse_cases().items() if c.get("class") in ("one_sided_reference", "upper_gate_3d")]
    rows = []
    for name in names:
        case = contract.parse_cases()[name]
        for mesh in (meshes or range(len(case["meshes"]))):
            problem, arrays = ref.case_problem(name, mesh)
            rows.append({"case": name, "mesh": mesh, "n_tets": int(len(arrays["tets"]))})
    Path(out).write_text(json.dumps(rows, indent=1) + "\n", encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("job", "worker"):
        p = sub.add_parser(name)
        p.add_argument("--case", required=True)
        p.add_argument("--mesh", type=int, required=True)
        p.add_argument("--formulation", required=True)
        p.add_argument("--solver", required=True, choices=SOLVERS)
        p.add_argument("--data", required=True)
        p.add_argument("--time-limit", type=float, default=TIME_LIMIT_S)
        p.add_argument("--mem-gb", type=float, default=MEMORY_GUARD_GB)
        p.add_argument("--overrides", default="{}")
        p.add_argument("--tag", default="")
        p.add_argument("--solve-only", action="store_true")
    p = sub.add_parser("queue")
    p.add_argument("--file", required=True)
    p.add_argument("--data", required=True)
    p.add_argument("--time-limit", type=float, default=TIME_LIMIT_S)
    p.add_argument("--mem-gb", type=float, default=MEMORY_GUARD_GB)
    p = sub.add_parser("table")
    p.add_argument("--data", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--md")
    args = parser.parse_args(argv)
    if args.command == "worker":
        worker(args.case, args.mesh, args.formulation, args.solver, args.data, args.time_limit, json.loads(args.overrides), args.tag, args.solve_only)
    elif args.command == "job":
        supervise(args.case, args.mesh, args.formulation, args.solver, args.data, args.time_limit, args.mem_gb, json.loads(args.overrides), args.tag, args.solve_only)
    elif args.command == "queue":
        queue(args.file, args.data, args.time_limit, args.mem_gb)
    elif args.command == "table":
        table(args.data, args.out, args.md)


if __name__ == "__main__":
    main()
