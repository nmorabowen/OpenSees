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
    python solver_study.py table  --data DIR --out SUMMARY.json [--md SUMMARY.md] [--base TEXT] [--diag DIR --diag-out FILE]
    python solver_study.py recheck --old DIR --new DIR [--only TEXT]  stored results through the current checker, no solve
    python solver_study.py sizes  --out FILE                           program sizes without solving (assemble only)
    python solver_study.py rankstudy --cases C:M ... --out FILE        dependent Lower equality rows, before and after presolve

Solvers: ``clarabel-qdldl``, ``clarabel-faer``, ``scs``, ``qoco`` (0.3.2; reported ``not_run`` when absent).
Presets (queue tag): see PRESETS; ``driver_presolve`` reduces and scales the program before the solver.
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

# Option variants named in queue files. ``driver_*`` keys are consumed by the driver, the rest go to Clarabel.
# presolve: singleton-row elimination, empty and proportional duplicate equality rows dropped, Ruiz scaling with
# cone-preserving row factors (see presolve_program). Calibration note in solver_study_results/presolve.md.
PRESETS = {"nodyn": {"dynamic_regularization_enable": False},
           "presolve": {"driver_presolve": True},
           "presolve_nodyn": {"driver_presolve": True, "dynamic_regularization_enable": False},
           "relgap1e-3": {"reduced_tol_gap_rel": 1e-3, "reduced_tol_gap_abs": 1e-3},
           "presolve_relgap1e-3": {"driver_presolve": True, "reduced_tol_gap_rel": 1e-3, "reduced_tol_gap_abs": 1e-3},
           "nodyn_relgap1e-3": {"dynamic_regularization_enable": False, "reduced_tol_gap_rel": 1e-3,
                                "reduced_tol_gap_abs": 1e-3}}
# Calibration of relgap1e-3 (2026-11-19, solver_study_results/presolve.md): Clarabel reports AlmostSolved when its
# last good iterate meets the reduced tolerances; the default reduced gap 5e-5 is missed by Lower runs that break
# down at gaps 5e-5 to 1e-3 with primal and dual residuals <= 2e-6. A looser gap changes only how close the stored
# stress field is to optimal (the Lower is lower by at most about the gap, 0.1 %); admissibility is proved by the
# checker on the full program, never by the solver status. It is not tuned to any reference value.
# Diagnostic variants: each changes one existing Clarabel option, for the cohesion Prandtl Lower stall study
# (2026-11-19). They are not study presets; only a variant that earns its place moves into PRESETS with a note.
DIAG_PRESETS = {
    "sreg1e-7": {"static_regularization_constant": 1e-7},
    "sreg1e-6": {"static_regularization_constant": 1e-6},
    "noequil": {"equilibrate_enable": False},
    "equil50": {"equilibrate_max_iter": 50, "equilibrate_min_scaling": 1e-6, "equilibrate_max_scaling": 1e6},
    "presolveoff": {"presolve_enable": False},
    "step95": {"max_step_fraction": 0.95},
    "irtight": {"iterative_refinement_reltol": 1e-15, "iterative_refinement_abstol": 1e-15,
                "iterative_refinement_max_iter": 50},
    "dyn1e-6": {"dynamic_regularization_delta": 1e-6, "dynamic_regularization_eps": 1e-10},
    "nodyn_sreg1e-7": {"dynamic_regularization_enable": False, "static_regularization_constant": 1e-7},
}
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


def solve_qoco(program, time_limit, overrides=None):
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
    # qoco 0.3.2: setup(n, m, p, P, c, A, b, G, h, l, nsoc, q, **settings); P=None is a linear objective. QOCO has
    # no time limit of its own; the supervisor's wall-clock kill is the budget.
    known = set(dir(qoco.qoco_ext.QOCOSettings))
    settings = {"verbose": 0}
    settings.update({k: v for k, v in (overrides or {}).items() if k in known})
    start = time.perf_counter()
    model.setup(n, G.shape[0], Aeq.shape[0], None, np.asarray(program.c, float), sp.csc_matrix(Aeq),
                np.asarray(program.b[:nz], float), sp.csc_matrix(G), np.asarray(program.b[nz:], float), 0, len(q),
                np.asarray(q, dtype=np.int32), **settings)
    solution = model.solve()
    wall = time.perf_counter() - start
    raw = str(solution.status)
    status = {"QOCO_SOLVED": "solved", "QOCO_SOLVED_INACCURATE": "almost_solved",
              "QOCO_MAX_ITER": "iteration_limit"}.get(raw, "numerical_error")
    out = {"status_raw": raw, "iterations": int(solution.iters), "solve_time_s": float(solution.solve_time_sec),
           "setup_time_s": float(solution.setup_time_sec), "obj_primal": float(solution.obj),
           "r_prim": float(solution.pres), "r_dual": float(solution.dres), "rel_gap": float(solution.gap)}
    used = {k: getattr(model.settings, k) for k in sorted(known) if not k.startswith("_")}
    from importlib.metadata import version
    record = {"name": "qoco", "version": version("qoco"), "settings": dict(used, status_raw=raw)}
    return status, np.asarray(solution.x, dtype=float), wall, record, out


DRIVER_KEYS = ("driver_presolve",)   # preset keys the driver consumes; never passed to a solver


def solve_program(program, solver, time_limit, overrides=None):
    overrides = dict(overrides or {})
    driver = {key: overrides.pop(key) for key in DRIVER_KEYS if key in overrides}
    if driver.get("driver_presolve"):
        start = time.perf_counter()
        reduced, stats = presolve_program(program)
        stats["wall_s"] = round(time.perf_counter() - start, 3)
        if reduced is None:     # the presolve refused: the solver gets the full program, and the row says so
            status, x, wall, record, info = solve_program(program, solver, time_limit, overrides)
        else:
            status, x, wall, record, info = solve_program(reduced, solver, time_limit, overrides)
            x = None if x is None else reduced.expand(x)
        info["driver_presolve"] = stats
        if record is not None:
            record["settings"]["driver_presolve"] = {"applied": reduced is not None, "fix_tol": PRESOLVE_FIX_TOL,
                                                     "parallel_tol": PRESOLVE_PARALLEL_TOL,
                                                     "ruiz_iters": PRESOLVE_RUIZ_ITERS}
        return status, x, wall, record, info
    if solver == "clarabel-qdldl":
        return solve_clarabel(program, "qdldl", time_limit, overrides)
    if solver == "clarabel-faer":
        return solve_clarabel(program, "faer", time_limit, overrides)
    if solver == "scs":
        return solve_scs(program, time_limit)
    if solver == "qoco":
        return solve_qoco(program, time_limit, overrides)
    raise SystemExit(f"unknown solver {solver}")


# ---------------------------------------------------------------------------------------------------
# Driver-side presolve: dependent equality rows and scaling (orchestrator request, 2026-11-19)
# ---------------------------------------------------------------------------------------------------
#
# The solver receives a reduced, scaled program; the field written for the checker is always the full-variable
# vector of the frozen program, so the checker still projects onto and verifies against every original row.
#
# Reduction, exact on the structure:
#   1. Singleton equality rows a_ij x_j = b_i fix x_j = b_i / a_ij. The fixed column is substituted into every row
#      (b <- b - A[:, j] x_j) and removed. Repeated until no singleton is left. Two singletons on one column must
#      agree to PRESOLVE_FIX_TOL (relative), else the presolve stops and the solver gets the full program.
#   2. Equality rows left with no entry are dependent; their right-hand side must be zero to PRESOLVE_FIX_TOL
#      (relative to the row's original norm and b), else the presolve stops as above.
#   3. Equality rows that are scalar multiples of an earlier row (same pattern, values proportional to the
#      relative tolerance) are dependent and dropped after the same right-hand-side check.
# Cone rows are never removed (a cone keeps its dimension); a fixed column only shifts their right-hand side.
# Any other linear dependence is left in place and is counted by ``rankstudy``.
#
# Scaling: PRESOLVE_RUIZ_ITERS sweeps of Ruiz equilibration in the infinity norm. Equality rows get one factor
# each; the rows of one second-order cone share one factor so that the cone is preserved; columns get one factor
# each. The solver sees A~ = R A C, b~ = R b, c~ = C c; the field is x = C x~. The objective value is unchanged.

PRESOLVE_FIX_TOL = 1e-12       # relative consistency of fixed values and of dropped right-hand sides
PRESOLVE_PARALLEL_TOL = 1e-12  # relative tolerance for proportional duplicate rows
PRESOLVE_RUIZ_ITERS = 10       # sweeps; calibration: row and column inf-norms within 2x of 1 after <= 10 sweeps


class Reduced:
    """A reduced program in the interface 15.4 shape, plus the map back to the full variable vector."""

    def __init__(self, A, b, c, n_zero, cones, free, fixed_values, n_full, col_scale, stats):
        self.A, self.b, self.c, self.n_zero, self.cones = A, b, c, n_zero, cones
        self.free, self.fixed_values, self.n_full, self.col_scale, self.stats = free, fixed_values, n_full, col_scale, stats

    def expand(self, x):
        import numpy as np
        full = self.fixed_values.copy()
        full[self.free] = self.col_scale * np.asarray(x, dtype=float)
        return full


def presolve_program(program, scale=True):
    import numpy as np
    import scipy.sparse as sp
    A = sp.csr_matrix(program.A, dtype=float)
    A.sum_duplicates()
    A.eliminate_zeros()
    b = np.asarray(program.b, dtype=float).copy()
    c = np.asarray(program.c, dtype=float).copy()
    m, n = A.shape
    nz = int(program.n_zero)
    row_norm0 = np.sqrt(np.asarray(A.multiply(A).sum(axis=1)).ravel())
    stats = {"rows": m, "n_zero": nz, "cols": n, "fixed_columns": 0, "empty_rows_dropped": 0,
             "duplicate_rows_dropped": 0, "sweeps": 0, "aborted": None}
    fixed = np.full(n, np.nan)
    alive_col = np.ones(n, dtype=bool)
    alive_row = np.ones(nz, dtype=bool)
    csc = A.tocsc()
    E = A[:nz].tocsr()
    counts = np.diff(E.indptr).astype(np.int64)          # live entries per equality row
    bref = np.maximum(row_norm0[:nz], np.abs(b[:nz]))
    # Singleton elimination with a work queue
    queue = [i for i in range(nz) if counts[i] == 1]
    Eind, Edat = E.indices, E.data
    while queue:
        stats["sweeps"] += 1
        next_queue = []
        for i in queue:
            if not alive_row[i] or counts[i] != 1:
                continue
            lo, hi = E.indptr[i], E.indptr[i + 1]
            live = [(Eind[k], Edat[k]) for k in range(lo, hi) if alive_col[Eind[k]]]
            if len(live) != 1:
                continue
            j, a = live[0]
            value = b[i] / a
            alive_row[i] = False
            if not alive_col[j]:
                continue
            fixed[j] = value
            alive_col[j] = False
            stats["fixed_columns"] += 1
            clo, chi = csc.indptr[j], csc.indptr[j + 1]
            rows, vals = csc.indices[clo:chi], csc.data[clo:chi]
            b[rows] -= vals * value
            for r in rows[rows < nz]:
                if alive_row[r]:
                    counts[r] -= 1
                    if counts[r] == 1:
                        next_queue.append(r)
                    elif counts[r] == 0:
                        alive_row[r] = False
                        if abs(b[r]) > PRESOLVE_FIX_TOL * max(1.0, bref[r]):
                            stats["aborted"] = f"inconsistent empty row {int(r)}: b={b[r]!r}"
                            return None, stats
                        stats["empty_rows_dropped"] += 1
        queue = next_queue
    c_fixed = np.where(alive_col, 0.0, np.nan_to_num(fixed))
    # Consistency of every removed equality row against the fixed values, on the original data.
    Afix = A[:nz][:, ~alive_col]
    removed = ~alive_row
    resid = (np.asarray(program.b, float)[:nz] - Afix @ fixed[~alive_col])[removed]
    # Rows removed may still have free columns only if they were singleton rows whose column got fixed elsewhere;
    # those rows have no free entries by construction.
    rel = np.abs(resid) / np.maximum(1.0, bref[removed])
    stats["removed_row_residual_max"] = float(rel.max()) if rel.size else 0.0
    if rel.size and rel.max() > PRESOLVE_FIX_TOL * 1e3:
        stats["aborted"] = f"removed rows inconsistent, max relative residual {rel.max()!r}"
        return None, stats
    free = np.flatnonzero(alive_col)
    keep_eq = np.flatnonzero(alive_row)
    Ared = A[:, free].tocsr()
    # Proportional duplicates among the surviving equality rows
    Eq = Ared[keep_eq]
    seen, drop = {}, []
    for local in range(Eq.shape[0]):
        lo, hi = Eq.indptr[local], Eq.indptr[local + 1]
        idx, val = Eq.indices[lo:hi], Eq.data[lo:hi]
        if idx.size == 0:
            continue
        order = np.argsort(idx)
        idx, val = idx[order], val[order]
        pivot = val[np.argmax(np.abs(val))]
        key = (tuple(idx.tolist()), tuple(np.round(val / pivot, 10).tolist()))
        row = keep_eq[local]
        if key in seen:
            first, first_pivot, first_val = seen[key]
            ratio = pivot / first_pivot
            if np.allclose(val, ratio * first_val, rtol=PRESOLVE_PARALLEL_TOL * 1e3, atol=0.0):
                if abs(b[row] - ratio * b[first]) > PRESOLVE_FIX_TOL * 1e3 * max(1.0, bref[row]):
                    stats["aborted"] = f"inconsistent duplicate rows {int(first)}, {int(row)}"
                    return None, stats
                drop.append(local)
                continue
        seen[key] = (row, pivot, val)
    stats["duplicate_rows_dropped"] = len(drop)
    keep_mask = np.ones(keep_eq.size, dtype=bool)
    keep_mask[drop] = False
    keep_eq = keep_eq[keep_mask]
    rows = np.concatenate([keep_eq, np.arange(nz, m)])
    Ared = Ared[rows].tocsr()
    bred = b[rows]
    cred = c[free]
    nz_red = int(keep_eq.size)
    # Ruiz scaling, cone-preserving
    col_scale = np.ones(free.size)
    if scale:
        row_scale = np.ones(Ared.shape[0])
        blocks = []
        start = nz_red
        for cone in program.cones:
            for _ in range(int(cone["count"])):
                blocks.append((start, start + int(cone["dim"])))
                start += int(cone["dim"])
        block_id = np.full(Ared.shape[0], -1)
        for k, (lo, hi) in enumerate(blocks):
            block_id[lo:hi] = k
        S = Ared.copy()
        for _ in range(PRESOLVE_RUIZ_ITERS):
            r = np.sqrt(np.asarray(abs(S).max(axis=1).todense()).ravel())
            if blocks:
                cone_rows = block_id >= 0
                bmax = np.zeros(len(blocks))
                np.maximum.at(bmax, block_id[cone_rows], r[cone_rows])
                r[cone_rows] = bmax[block_id[cone_rows]]
            r[r == 0.0] = 1.0
            cmax = np.sqrt(np.asarray(abs(S).max(axis=0).todense()).ravel())
            cmax[cmax == 0.0] = 1.0
            S = sp.diags(1.0 / r) @ S @ sp.diags(1.0 / cmax)
            row_scale /= r
            col_scale /= cmax
        Ared = S.tocsr()
        bred = row_scale * bred
        cred = col_scale * cred
        stats["row_scale_range"] = [float(row_scale.min()), float(row_scale.max())]
        stats["col_scale_range"] = [float(col_scale.min()), float(col_scale.max())]
    norms = np.sqrt(np.asarray(Ared[:nz_red].multiply(Ared[:nz_red]).sum(axis=1)).ravel())
    stats.update(rows_out=int(Ared.shape[0]), n_zero_out=nz_red, cols_out=int(free.size),
                 eq_rows_removed=int(nz - nz_red), objective_offset=float(np.asarray(program.c, float) @ c_fixed),
                 eq_row_norm_spread_in=float(row_norm0[:nz][row_norm0[:nz] > 0].max() / row_norm0[:nz][row_norm0[:nz] > 0].min()),
                 eq_row_norm_spread_out=float(norms.max() / norms.min()) if norms.size else None)
    fixed_values = np.nan_to_num(fixed)
    fixed_values[free] = 0.0
    return Reduced(Ared, bred, cred, nz_red, list(program.cones), free, fixed_values, n, col_scale, stats), stats


def equality_nullity(A_eq, shifts=(1e-10, 1e-9)):
    """Number of dependent rows of ``A_eq``: the count of eigenvalues of M = A A^T below a small shift.

    M - sigma I is factorised by SuperLU with no row pivoting in symmetric mode; by Sylvester's law of inertia the
    number of negative pivots is the number of eigenvalues below sigma. Rows with no entry are counted directly.
    The count is taken at two shifts (relative to the largest row norm squared, an upper bound on nothing but a
    scale); equal counts mean a clear spectral gap. Returns (count at each shift, rows with entries, empty rows).
    """
    import numpy as np
    import scipy.sparse as sp
    import scipy.sparse.linalg as sla
    A = sp.csr_matrix(A_eq)
    nnz_rows = np.diff(A.indptr) > 0
    A = A[nnz_rows]
    M = (A @ A.T).tocsc()
    scale = float(abs(M).max())
    counts = []
    for shift in shifts:
        sigma = shift * scale
        lu = sla.splu(M - sigma * sp.identity(M.shape[0], format="csc"), permc_spec="MMD_AT_PLUS_A",
                      diag_pivot_thresh=0.0, options={"SymmetricMode": True})
        pivots = lu.U.diagonal()
        if not np.array_equal(lu.perm_r, lu.perm_c):
            counts.append(None)          # row pivoting happened: the inertia count does not hold
            continue
        counts.append(int(np.sum(pivots < 0.0)))
    return counts, int(nnz_rows.sum()), int((~nnz_rows).sum())


def rank_study(cases, out):
    """Per case and mesh: equality rows, dependent rows before and after the presolve, row-norm spread."""
    import numpy as np
    import limit_reference as ref
    rows = []
    for spec in cases:
        name, mesh = spec.rsplit(":", 1)
        mesh = int(mesh)
        start = time.perf_counter()
        problem, arrays = ref.case_problem(name, mesh)
        program = ref.assemble(problem, arrays, "lower.p1d-vertex.v0", "q")
        nz = int(program.n_zero)
        import scipy.sparse as sp
        A = sp.csr_matrix(program.A)
        norms = np.sqrt(np.asarray(A[:nz].multiply(A[:nz]).sum(axis=1)).ravel())
        item = {"case": name, "mesh": mesh, "n_tets": int(len(arrays["tets"])), "eq_rows": nz, "cols": int(A.shape[1]),
                "eq_row_norm_min": float(norms[norms > 0].min()), "eq_row_norm_max": float(norms.max())}
        try:
            counts, live, empty = equality_nullity(A[:nz])
            item.update(dependent_before=counts, empty_rows_before=empty)
        except Exception as error:
            item["dependent_before"] = f"{type(error).__name__}: {error}"[:200]
        reduced, stats = presolve_program(program, scale=False)
        item["presolve"] = stats
        if reduced is not None:
            try:
                counts, live, empty = equality_nullity(reduced.A[:reduced.n_zero])
                item.update(dependent_after=counts, empty_rows_after=empty)
            except Exception as error:
                item["dependent_after"] = f"{type(error).__name__}: {error}"[:200]
        item["wall_s"] = round(time.perf_counter() - start, 2)
        item["peak_gb"] = round(peak_gb(), 3)
        rows.append(item)
        print(json.dumps({k: v for k, v in item.items() if k != "presolve"}), json.dumps(stats), flush=True)
        Path(out).write_text(json.dumps(rows, indent=1) + "\n", encoding="utf-8")
    return rows


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
    solver_overrides = {k: v for k, v in (overrides or {}).items() if k not in DRIVER_KEYS}
    if solver_overrides:
        patch_interior_solves(solver_overrides)
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
        # One file name per solver so that solvers on one problem do not overwrite each other.
        directory = data / "certs" / (solver + (f".{tag}" if tag else ""))
        directory.mkdir(parents=True, exist_ok=True)
        # The result files are large; keep them only in the certificate directory, beside the gate record.
        (data / f"{stem_ref}.result.json").replace(directory / f"{stem_ref}.result.json")
        (data / f"{stem_ref}.result.npz").replace(directory / f"{stem_ref}.result.npz")
        # The checker reads the files the reference wrote, as an independent consumer would.
        certify(row, data / f"{ref.problem_stem(problem)}.problem.json", directory / f"{stem_ref}.result.json",
                directory, data)
        phase("check", t)
    except MemoryError as error:
        row["error"] = {"type": "MemoryError", "message": str(error)[:300], "trace": traceback.format_exc()[-1500:]}
    except Exception as error:
        row["error"] = {"type": type(error).__name__, "message": str(error)[:300], "trace": traceback.format_exc()[-1500:]}
    finally:
        row["peak_gb"]["final"] = round(peak_gb(), 3)
        save()


def certify(row, problem_path, result_path, directory, data):
    """Run the independent checker on one stored problem and result and write what it returns.

    ``write_outcome`` writes ``<stem>.certificate.json`` with its ``.certificate.method.json`` beside it, or a
    ``<stem>.refusal.json`` instead (interface 15.6, 15.8, 15.9). Gate admissibility is decided only by
    ``read_gate_certificate`` on the written file; the method report is also checked once more against the result
    arrays, so that its interior-field hash is bound to the field actually stored.
    """
    import limit_certificate as lc
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for stale in ("certificate_file", "gate_admissible", "gate_refusal", "interior_field", "kappa_method",
                  "kappa_eq", "repair_mu", "method_failures_with_arrays"):
        row.pop(stale, None)
    cproblem, carrays = lc.load_document(problem_path)
    cresult, cfields = lc.load_document(result_path)
    try:
        outcome = lc.check(cproblem, carrays, cresult, cfields)
    except lc.Refusal as error:      # malformed documents: nothing is written
        row["certificate"] = {"classification": "refusal", "code": error.code, "note": f"{error}"[:400]}
        return None
    path = lc.write_outcome(directory, cresult, outcome)
    row["certificate_file"] = str(path.relative_to(data))
    if outcome.refusal is not None:
        refusal = outcome.refusal
        row["certificate"] = {"classification": "capacity_refusal", "kind": refusal["kind"], "stage": refusal["stage"],
                              "note": refusal["detail"][:400]}
        row["gate_admissible"] = False
        return path
    cert = outcome.certificate
    row["certificate"] = {
        "classification": cert["classification"], "reasons": list(cert["reasons"]),
        "value_candidate": cert["value_candidate"], "value_certified": cert["value_certified"], "floor": cert["floor"],
        "domain_claim": cert["domain_claim"], "criteria": cert["criteria"], "flow_rule": cert["flow_rule"],
        "interface_laws": cert["interface_laws"], "checks": cert["checks"],
        "method_report": outcome.method is not None,
    }
    row["kappa_eq"] = cert["checks"]["kappa_eq"]
    row["repair_mu"] = cert["checks"]["repair_mu"]
    if outcome.method is not None:
        row["interior_field"] = outcome.method["interior_field"]
        row["kappa_method"] = outcome.method["kappa"]["method"]
        row["method_failures_with_arrays"] = lc.method_failures(cert, outcome.method, cfields)
    try:
        lc.read_gate_certificate(path)
        row["gate_admissible"] = not row.get("method_failures_with_arrays")
        if not row["gate_admissible"]:
            row["gate_refusal"] = row["method_failures_with_arrays"][0]
    except lc.Refusal as error:
        row["gate_admissible"] = False
        row["gate_refusal"] = str(error)[:300]
    return path


def normalise_outcome(row):
    """A job counts as completed only when the supervisor saw a clean exit and the worker raised nothing.

    A worker exception (for example ``ModuleNotFoundError``) is caught and stored in ``error`` with exit code 0;
    that row is an ``error``, never ``completed``. A row that the supervisor never finished is ``incomplete``.
    """
    if row.get("error"):
        row["outcome"] = "error"
        row.setdefault("certificate", {"classification": "no_value", "note": f"error {row['error'].get('type')}"})
    elif row.get("outcome") is None:
        row["outcome"] = "incomplete"
        row.setdefault("certificate", {"classification": "no_value", "note": "supervisor record missing"})
    return row


FAILED_OUTCOMES = ("error", "crash", "incomplete", "memory_guard", "timeout")


# ---------------------------------------------------------------------------------------------------
# Supervisor
# ---------------------------------------------------------------------------------------------------

def run_guarded(command, wall_limit, mem_gb):
    """Run one child process under the wall-clock limit and the resident-memory guard."""
    import psutil
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
        if time.perf_counter() - start > wall_limit:
            outcome = "timeout"
            process.kill()
            break
        time.sleep(POLL_S)
    try:
        tail = process.communicate(timeout=30)[0][-600:]
    except Exception:
        tail = ""
    return outcome, time.perf_counter() - start, peak, process.returncode, tail


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
    outcome, wall, peak, returncode, tail = run_guarded(command, time_limit + SUPERVISOR_SLACK_S, mem_gb)
    row = json.loads(row_path.read_text(encoding="utf-8")) if row_path.exists() else {
        "format": ROW_FORMAT, "format_version": ROW_VERSION, "case": case, "mesh": mesh,
        "formulation": formulation, "solver": solver, "phases": {}, "peak_gb": {}}
    row.update(outcome=outcome, supervisor_wall_s=round(wall, 2), supervisor_peak_rss_gb=round(peak, 3),
               exit_code=returncode, memory_guard_gb=mem_gb)
    if outcome == "completed" and returncode not in (0, None):
        row["outcome"] = "crash"
        row["crash_tail"] = tail
    if outcome != "completed":
        row.setdefault("certificate", {"classification": "no_value", "note": outcome})
    normalise_outcome(row)
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
            supervise(line[0], int(line[1]), line[2], line[3], data, time_limit, mem_gb, PRESETS.get(tag, DIAG_PRESETS.get(tag)), tag)


# ---------------------------------------------------------------------------------------------------
# Recheck: stored results under an amended checker, without solving again
# ---------------------------------------------------------------------------------------------------

RECHECK_TIME_LIMIT_S = 1800.0   # one checker run; the slowest stored check took about 50 s


def stored_result(old, row):
    """The stored result of one row in the old data directory, or None.

    A certified or estimated row names its certificate; the result sits beside it. A row stopped after the result
    was written and before the checker ran leaves the result in the top directory; it is taken only when its solver
    record names the row's solver path and option overrides.
    """
    file = row.get("certificate_file")
    if file:
        path = old / Path(file.replace("\\", "/"))
        result = path.with_name(path.name[: -len(".certificate.json")] + ".result.json")
        return result if result.exists() else None
    stem = f"{row['case']}.m{row['mesh']}.q.{row['formulation']}"
    path = old / f"{stem}.result.json"
    if not path.exists():
        return None
    record = json.loads(path.read_text(encoding="utf-8")).get("solver", {})
    settings = record.get("settings", {})
    if row["solver"].startswith("clarabel"):
        same = record.get("name") == "clarabel" and settings.get("direct_solve_method") == row["solver"].split("-")[1] \
            and settings.get("overrides", {}) == row.get("overrides", {})
    else:
        same = record.get("name") == row["solver"]
    return path if same else None


def recheck_one(old, new, row_name):
    """Worker of ``recheck``: the amended checker on one stored problem and result. Writes the new row."""
    old, new = Path(old), Path(new)
    row = normalise_outcome(json.loads((old / row_name).read_text(encoding="utf-8")))
    previous = row.get("certificate", {})
    row["old_certificate"] = {"classification": previous.get("classification"),
                              "value_certified": previous.get("value_certified"),
                              "reasons": previous.get("reasons"), "file": row.get("certificate_file")}
    result = stored_result(old, row)
    row["recheck"] = {"result_file": None if result is None else str(result)}
    start = time.perf_counter()
    try:
        if result is None:
            row.pop("certificate_file", None)
            row["certificate"] = dict(previous) if previous.get("classification") == "no_value" else \
                {"classification": "no_value", "note": "no stored result"}
        else:
            label = row["solver"] + (f".{row['tag']}" if row.get("tag") else "")
            problem = old / f"{row['case']}.m{row['mesh']}.problem.json"
            certify(row, problem, result, new / "certs" / label, new)
    except Exception as error:
        row["recheck"]["error"] = {"type": type(error).__name__, "message": str(error)[:300],
                                   "trace": traceback.format_exc()[-1500:]}
        row["certificate"] = {"classification": "no_check", "note": f"recheck error {type(error).__name__}"}
    row["recheck"]["wall_s"] = round(time.perf_counter() - start, 3)
    row["recheck"]["peak_gb"] = round(peak_gb(), 3)
    (new / row_name).write_text(json.dumps(row, indent=1, sort_keys=True) + "\n", encoding="utf-8")


def recheck(old, new, mem_gb=MEMORY_GUARD_GB, only=None):
    """Every stored row of ``old`` through the current checker into ``new``; ``old`` is only read."""
    old, new = Path(old), Path(new)
    new.mkdir(parents=True, exist_ok=True)
    names = sorted(p.name for p in old.glob("*.row.json"))
    for name in names:
        if only and only not in name:
            continue
        target = new / name
        if target.exists():
            target.unlink()
        command = [sys.executable, str(Path(__file__).resolve()), "recheck-one", "--old", str(old), "--new", str(new),
                   "--row", name]
        outcome, wall, peak, code, tail = run_guarded(command, RECHECK_TIME_LIMIT_S, mem_gb)
        if target.exists():
            row = json.loads(target.read_text(encoding="utf-8"))
        else:   # the guard stopped the checker: the old row, with no certificate
            row = normalise_outcome(json.loads((old / name).read_text(encoding="utf-8")))
            previous = row.get("certificate", {})
            row["old_certificate"] = {"classification": previous.get("classification"),
                                      "value_certified": previous.get("value_certified"),
                                      "reasons": previous.get("reasons"), "file": row.pop("certificate_file", None)}
            row["certificate"] = {"classification": "no_check", "note": f"recheck {outcome} exit {code}"}
            row["recheck"] = {"tail": tail}
        row["recheck"].update(outcome=outcome if code in (0, None) or outcome != "completed" else "crash",
                              supervisor_wall_s=round(wall, 2), supervisor_peak_rss_gb=round(peak, 3), exit_code=code)
        target.write_text(json.dumps(row, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        old_c, new_c = row["old_certificate"]["classification"], row["certificate"].get("classification")
        print(f"{name}: {row['recheck']['outcome']} {old_c} -> {new_c} wall={wall:.1f}s peak={peak:.2f}GB "
              f"{row.get('gate_refusal', '')}", flush=True)


# ---------------------------------------------------------------------------------------------------
# Table and sizes
# ---------------------------------------------------------------------------------------------------

def _midpoint_gap(low, high):
    """The skill's gap: (U - L) / ((U + L) / 2), only for a certified pair with 0 < L <= U."""
    return (high - low) / (0.5 * (high + low)) if low is not None and high is not None and 0.0 < low <= high else None


def gate_records(data, rows, case):
    """The checker records of one case for the frozen contract evaluator (interface 17.9).

    Per (mesh, formulation) the best gate-admissible certificate is offered (largest Lower, smallest Upper), with its
    result document and method report; where none is admissible, every estimate certificate and refusal is offered
    so that the evaluator names the cause. Result documents are read from the stored result files.
    """
    import limit_certificate as lc
    chosen, others = {}, {}
    for row in rows:
        if row["case"] != case or not row.get("certificate_file") or not row.get("recheck", {}).get("result_file"):
            continue
        key = (row["mesh"], row["formulation"])
        if row.get("gate_value") is not None:
            current = chosen.get(key)
            value = row["gate_value"]
            if current is None or (value > current["gate_value"] if row["bound"] == "lower" else value < current["gate_value"]):
                chosen[key] = row
        else:
            others.setdefault(key, []).append(row)
    offered = list(chosen.values()) + [row for key, items in others.items() if key not in chosen for row in items]
    certificates, results, methods, refusals = [], [], [], []
    for row in offered:
        path = data / row["certificate_file"]
        results.append(json.loads(Path(row["recheck"]["result_file"]).read_text(encoding="utf-8")))
        if path.name.endswith(".refusal.json"):
            refusals.append(lc.load_document(path)[0])
            continue
        certificate, method = lc.read_certificate(path)
        certificates.append(certificate)
        if method is not None:
            methods.append(method)
    return certificates, results, methods, refusals, chosen


def gate_verdicts(data, rows):
    """Every active case with a frozen gap gate, judged by the contract's own evaluator on the rechecked records."""
    import test_limit_contract as contract
    out = []
    for name in contract.ACTIVE:
        case = contract.ALL_CASES[name]
        gap_gates = [tokens for tokens in case["gate"] if tokens[0] == "gap_le"]
        if not gap_gates:
            continue
        limit_name, gate_mesh = gap_gates[0][1], int(gap_gates[0][3])
        certificates, results, methods, refusals, chosen = gate_records(data, rows, name)
        verdict = contract.evaluate_case(case, certificates, results, (), methods, refusals=refusals)
        lowers = [f for f in case["formulation"] if f.startswith("lower")]
        uppers = [f for f in case["formulation"] if f.startswith("upper")] + list(case.get("control") or [])
        mine = [row for row in rows if row["case"] == name]
        meshes_run = sorted({row["mesh"] for row in mine})
        per_mesh = []
        for mesh in range(len(case["mesh"])):
            entry = {"mesh": mesh, "run": mesh in meshes_run}
            low = chosen.get((mesh, lowers[0]))
            entry["lower"] = None if low is None else {"value": low["gate_value"], "run": low["label"]}
            if low is None:
                estimates = [r["value_estimate"] for r in mine if r["mesh"] == mesh and r["formulation"] == lowers[0]
                             and r.get("value_estimate") is not None]
                if estimates:
                    entry["lower_estimate"] = max(estimates)
            for upper in uppers:
                kind = "jump" if "jump" in upper else "control"
                up = chosen.get((mesh, upper))
                entry[f"upper_{kind}"] = None if up is None else {"value": up["gate_value"], "run": up["label"]}
                entry[f"gap_{kind}"] = _midpoint_gap(None if low is None else low["gate_value"],
                                                     None if up is None else up["gate_value"])
                if up is None:
                    estimates = [r["value_estimate"] for r in mine if r["mesh"] == mesh and r["formulation"] == upper
                                 and r.get("value_estimate") is not None]
                    if estimates:
                        entry[f"upper_{kind}_estimate"] = min(estimates)
            per_mesh.append(entry)
        gap_items = [item for item in verdict["items"] if item["kind"] == "gap_le"]
        out.append({
            "case": name, "class": case["class"], "gap_limit": limit_name,
            "gap_limit_value": contract.TOLERANCES[limit_name], "gate_mesh": gate_mesh,
            "gate_mesh_run": gate_mesh in meshes_run, "meshes_run": meshes_run, "gap_gate_items": gap_items,
            "gap_gate_passed": bool(gap_items) and all(item["status"] == "pass" for item in gap_items),
            "case_passed": verdict["passed"], "claim_withheld": verdict["claim_withheld"],
            "p0_blocking": verdict["p0_blocking"], "stop_program": verdict["stop_program"],
            "invalid_records": verdict["invalid_records"],
            "failed_items": [[item["kind"], item["mesh"], item["formulation"], item["cause"]] for item in verdict["items"]
                             if item["status"] == "fail"],
            "per_mesh": per_mesh,
        })
    return out


def table(data, out, md=None, base=""):
    """Summary of every row, with the best certified bracket per (case, mesh) and the frozen gap gates.

    Certified values are read only through ``read_gate_certificate``. The best Lower is the largest certified Lower
    of any solver run, the best Upper the smallest certified Upper of one formulation; each is a certified bound on
    its own, so the pair is a bracket as the checker's ``bracket`` decides, never as arithmetic here. Gate verdicts
    come from the contract's own ``evaluate_case``. A row whose job raised, crashed, ran out of memory or time, or
    never finished counts as a failed job, whatever its certificate says.
    """
    import datetime
    import limit_certificate as lc
    data = Path(data)
    rows = [normalise_outcome(json.loads(p.read_text(encoding="utf-8"))) for p in sorted(data.glob("*.row.json"))]
    gate = {}
    for index, row in enumerate(rows):
        row["label"] = row["solver"] + (f".{row['tag']}" if row.get("tag") else "")
        file = row.get("certificate_file")
        row["gate_value"] = None
        if file and file.endswith(".certificate.json") and not row.get("method_failures_with_arrays"):
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
    changes = []
    for row in rows:
        old = row.get("old_certificate")
        if old is None:
            continue
        before, after = old.get("classification"), row.get("certificate", {}).get("classification")
        if before != after or (before == "certified_bound") != (row["gate_value"] is not None):
            changes.append({"row": stem_of(row["case"], row["mesh"], row["formulation"], row["solver"], row.get("tag", "")),
                            "old": before, "new": after, "gate_admissible": row["gate_value"] is not None,
                            "old_value": old.get("value_certified"), "new_value": row["gate_value"],
                            "reasons": row.get("certificate", {}).get("reasons"), "refusal": row.get("gate_refusal"),
                            "note": row.get("certificate", {}).get("note")})
    failures = [stem_of(r["case"], r["mesh"], r["formulation"], r["solver"], r.get("tag", "")) + f": {r['outcome']}"
                for r in rows if r["outcome"] in FAILED_OUTCOMES]
    recheck_failures = [stem_of(r["case"], r["mesh"], r["formulation"], r["solver"], r.get("tag", ""))
                        + f": recheck {r['recheck'].get('outcome')}" for r in rows
                        if r.get("recheck") and r["recheck"].get("outcome") not in (None, "completed")]
    gates = gate_verdicts(data, rows)
    summary = {"format": "ladruno-limit-solver-study-summary", "format_version": 1,
               "date": datetime.date.today().isoformat(), "reference_checker_base": base, "data": str(data),
               "driver": "Ladruno_files/testbed/limit_analysis/solver_study.py",
               "counts": {"rows": len(rows), "failed_jobs": len(failures), "failed_rechecks": len(recheck_failures),
                          "certified_gate_admissible": sum(r["gate_value"] is not None for r in rows),
                          "classification_changes": len(changes)},
               "failed_jobs": failures, "failed_rechecks": recheck_failures, "classification_changes": changes,
               "gates": gates, "rows": rows, "brackets": brackets}
    Path(out).write_text(json.dumps(summary, indent=1, sort_keys=True) + chr(10), encoding="utf-8")
    if md:
        def fmt(value):
            return "-" if value is None else f"{value:.6g}"
        lines = ["# WP-237d solver study, rechecked under the amended checker", "",
                 f"Base: {base}. Driver: `Ladruno_files/testbed/limit_analysis/solver_study.py` (`table`). "
                 f"Data: `{data}`. Date: {summary['date']}. Rows {len(rows)}; failed jobs {len(failures)}; failed "
                 f"rechecks {len(recheck_failures)}; gate-admissible certificates "
                 f"{summary['counts']['certified_gate_admissible']}; classification changes {len(changes)}.", "",
                 "Every value below is a certified bound read through `read_gate_certificate`, unless it is in the "
                 "estimates column. Gap = (U - L) / ((U + L) / 2), computed only for a certified pair on one mesh.", "",
                 "## Frozen gap gates (contract `evaluate_case`)", "",
                 "| case | gate | limit | gate mesh run | meshes run | gap gate | case passed | withheld | P0 blocking |",
                 "|---|---|---|---|---|---|---|---|---|"]
        for g in gates:
            causes = sorted({item["cause"] for item in g["gap_gate_items"] if item["status"] == "fail"})
            lines.append(f"| {g['case']} | {g['gap_limit']} m{g['gate_mesh']} | {g['gap_limit_value']} | "
                         f"{g['gate_mesh_run']} | {g['meshes_run']} | "
                         f"{'pass' if g['gap_gate_passed'] else 'fail: ' + ', '.join(causes)} | {g['case_passed']} | "
                         f"{g['claim_withheld']} | {g['p0_blocking']} |")
        lines += ["", "## Best certified bounds per mesh", "",
                  "| case | m | Lower | Upper jump | gap jump | Upper P2c | gap P2c | estimates (no certificate) |",
                  "|---|---|---|---|---|---|---|---|"]
        for g in gates:
            for e in g["per_mesh"]:
                if not e["run"]:
                    continue
                est = [f"{key.replace('_estimate', '')}={e[key]:.6g}" for key in
                       ("lower_estimate", "upper_jump_estimate", "upper_control_estimate") if key in e]
                lines.append(f"| {g['case']} | {e['mesh']} | {fmt((e['lower'] or {}).get('value'))} | "
                             f"{fmt((e.get('upper_jump') or {}).get('value'))} | {fmt(e.get('gap_jump'))} | "
                             f"{fmt((e.get('upper_control') or {}).get('value'))} | {fmt(e.get('gap_control'))} | "
                             f"{'; '.join(est)} |")
        lines += ["", "## Classification changes under the amended checker", ""]
        lines += [f"- `{c['row']}`: {c['old']} -> {c['new']}, gate admissible {c['gate_admissible']}, reasons "
                  f"{c['reasons']} {c['refusal'] or ''} {c['note'] or ''}".rstrip() for c in changes] or ["- none"]
        lines += ["", "## Failed jobs (error, crash, incomplete, memory guard, timeout)", ""]
        lines += [f"- `{f}`" for f in failures + recheck_failures] or ["- none"]
        lines += ["", "## Rows", "",
                  "| case | m | formulation | solver | outcome | status | it | solve s | peak GB | old cert | cert | value | kappa_eq |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in rows:
            info = r.get("solver_info", {})
            lines.append("| {case} | {mesh} | {f} | {solver} | {o} | {s} | {it} | {w} | {pk} | {oc} | {c} | {v} | {k} |".format(
                case=r["case"], mesh=r["mesh"], f=r["formulation"].split(".")[1], solver=r["label"],
                o=r.get("outcome"), s=r.get("solver_status"), it=info.get("iterations"),
                w=r.get("solver_wall_s"), pk=r.get("supervisor_peak_rss_gb"),
                oc=(r.get("old_certificate") or {}).get("classification"),
                c=r.get("certificate", {}).get("classification"),
                v=None if r.get("gate_value") is None else f"{r['gate_value']:.6g}", k=r.get("kappa_eq")))
        Path(md).write_text(chr(10).join(lines) + chr(10), encoding="utf-8")
    return summary


def diagnostics(data, out, base=""):
    """The solve-only diagnostic rows (no certificate is attempted), with job outcomes normalised."""
    import datetime
    rows = [normalise_outcome(json.loads(p.read_text(encoding="utf-8"))) for p in sorted(Path(data).glob("*.row.json"))]
    failures = [stem_of(r["case"], r["mesh"], r["formulation"], r["solver"], r.get("tag", "")) + f": {r['outcome']}"
                for r in rows if r["outcome"] in FAILED_OUTCOMES]
    doc = {"format": "ladruno-limit-solver-study-diagnostics", "date": datetime.date.today().isoformat(),
           "reference_checker_base": base, "data": str(data), "failed_jobs": failures, "rows": rows}
    Path(out).write_text(json.dumps(doc, indent=1, sort_keys=True) + chr(10), encoding="utf-8")


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
    p.add_argument("--base", default="", help="reference and checker commits the certificates were made under")
    p.add_argument("--diag", help="directory of solve-only diagnostic rows")
    p.add_argument("--diag-out")
    p = sub.add_parser("recheck")
    p.add_argument("--old", required=True, help="data directory holding the stored rows and results (read only)")
    p.add_argument("--new", required=True, help="new data directory for rows, certificates and method reports")
    p.add_argument("--mem-gb", type=float, default=MEMORY_GUARD_GB)
    p.add_argument("--only", help="substring of the row names to recheck")
    p = sub.add_parser("rankstudy")
    p.add_argument("--cases", nargs="+", required=True, help="case:mesh items")
    p.add_argument("--out", required=True)
    p = sub.add_parser("recheck-one")
    p.add_argument("--old", required=True)
    p.add_argument("--new", required=True)
    p.add_argument("--row", required=True)
    args = parser.parse_args(argv)
    if args.command == "worker":
        worker(args.case, args.mesh, args.formulation, args.solver, args.data, args.time_limit, json.loads(args.overrides), args.tag, args.solve_only)
    elif args.command == "job":
        supervise(args.case, args.mesh, args.formulation, args.solver, args.data, args.time_limit, args.mem_gb, json.loads(args.overrides), args.tag, args.solve_only)
    elif args.command == "queue":
        queue(args.file, args.data, args.time_limit, args.mem_gb)
    elif args.command == "table":
        table(args.data, args.out, args.md, args.base)
        if args.diag:
            diagnostics(args.diag, args.diag_out, args.base)
    elif args.command == "recheck":
        recheck(args.old, args.new, args.mem_gb, args.only)
    elif args.command == "rankstudy":
        rank_study(args.cases, args.out)
    elif args.command == "recheck-one":
        recheck_one(args.old, args.new, args.row)


if __name__ == "__main__":
    main()
