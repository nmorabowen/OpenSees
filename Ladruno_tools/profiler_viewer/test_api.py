"""
test_api.py  —  plain-assert guard for the FastAPI backend (profiler_api.py, P7).

Builds a small two-run fixture (one with a per-step series, one without) with the
same write_run() contract the C++ engine mirrors, then exercises every endpoint
through Starlette's TestClient (no live server needed). Run:

    python test_api.py
"""

from __future__ import annotations

import os
import tempfile

import h5py
from fastapi.testclient import TestClient

from profiler_schema import write_run
from make_sample import base_meta, implicit_rollup, make_series
from profiler_api import create_app


def _write_fixture(path: str) -> None:
    with h5py.File(path, "w") as f:
        write_run(f, "caseA",
                  meta=base_meta(integrator="Newmark", algorithm="Newton", solver="Mumps",
                                 nSteps=10, dt_min=0.005, dt_max=0.005, wall_ms_total=1000.0),
                  rollup=implicit_rollup(form_tangent_ms=400.0, factor_ms=250.0),
                  series=make_series(10, 0.005,
                                     {"formTangent": 3.3, "linearSolve": 2.2,
                                      "update": 0.8, "formUnbalance": 1.0, "convTest": 0.05}),
                  memory=dict(matrix_live=160, vector_live=384, id_live=732, peak_bytes=2844,
                              components_live=[(12, 20), (1, 72)]))
        # caseB: no series -> exercises the /series 404 path
        write_run(f, "caseB",
                  meta=base_meta(integrator="Newmark", algorithm="ModifiedNewton", solver="Mumps",
                                 nSteps=10, dt_min=0.005, dt_max=0.005, wall_ms_total=600.0),
                  rollup=implicit_rollup(form_tangent_ms=130.0, factor_ms=75.0),
                  series=None,
                  memory=dict(matrix_live=160, vector_live=384, id_live=732, peak_bytes=2844,
                              components_live=[(12, 20), (1, 72)]))


def main() -> int:
    problems = 0

    def check(cond, msg):
        nonlocal problems
        if cond:
            print(f"  [OK] {msg}")
        else:
            print(f"  [FAIL] {msg}")
            problems += 1

    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "profile_test.h5")
        _write_fixture(path)
        client = TestClient(create_app(path))

        # -- /health --------------------------------------------------------
        r = client.get("/health")
        check(r.status_code == 200, "/health 200")
        h = r.json()
        check(h["status"] == "ok", "/health status ok")
        check(set(h["runs"]) == {"caseA", "caseB"}, f"/health lists both runs (got {h['runs']})")

        # -- /runs (manifest) -----------------------------------------------
        r = client.get("/runs")
        check(r.status_code == 200, "/runs 200")
        rows = r.json()
        check(len(rows) == 2, f"/runs has 2 rows (got {len(rows)})")
        check({row["run"] for row in rows} == {"caseA", "caseB"}, "/runs covers both runs")
        check(all("wall_ms_total" in row for row in rows), "/runs rows carry wall_ms_total")

        # -- /runs/{run}/header ---------------------------------------------
        r = client.get("/runs/caseA/header")
        check(r.status_code == 200, "/runs/caseA/header 200")
        check(r.json()["integrator"] == "Newmark", "header integrator == Newmark")

        # -- /runs/{run}/rollup ---------------------------------------------
        r = client.get("/runs/caseA/rollup")
        check(r.status_code == 200, "/runs/caseA/rollup 200")
        root = r.json()
        check(root.get("children"), "rollup root has children")
        # the P0#1 per-element breakdown should survive the round trip
        scs = next((c for c in root["children"] if c["name"] == "solveCurrentStep"), None)
        ft = next((c for c in (scs["children"] if scs else []) if c["name"] == "formTangent"), None)
        check(ft is not None and ft.get("elem_by_type"),
              "rollup formTangent carries elem_by_type (P0#1)")

        # -- /runs/{run}/series ---------------------------------------------
        r = client.get("/runs/caseA/series")
        check(r.status_code == 200, "/runs/caseA/series 200")
        s = r.json()
        check(len(s["step"]) == 10 and "phases" in s, "series has 10 steps + phases")
        # caseB has no series -> 404
        r = client.get("/runs/caseB/series")
        check(r.status_code == 404, f"/runs/caseB/series 404 (no per-step series) (got {r.status_code})")

        # -- /runs/{run}/memory ---------------------------------------------
        r = client.get("/runs/caseA/memory")
        check(r.status_code == 200, "/runs/caseA/memory 200")
        m = r.json()
        check(m["id_live"] == 732, "memory id_live round-trips")
        tags = {c["classTag"] for c in m["components_live"]}
        check(tags == {12, 1}, f"memory census classTags round-trip (got {tags})")

        # -- /diff ----------------------------------------------------------
        r = client.get("/diff", params={"base": "caseA", "cand": "caseB"})
        check(r.status_code == 200, "/diff 200")
        diff = r.json()
        check(diff["base"] == "caseA" and diff["cand"] == "caseB", "/diff echoes base/cand")
        check(len(diff["rows"]) > 0, "/diff returns rows")
        ft_row = next((row for row in diff["rows"] if row["name"] == "formTangent"), None)
        check(ft_row is not None and ft_row["d_abs_ms"] < 0,
              "/diff shows formTangent faster in caseB (negative d_abs_ms)")

        # -- error paths ----------------------------------------------------
        check(client.get("/runs/nope/rollup").status_code == 404, "unknown run -> 404 (rollup)")
        check(client.get("/runs/nope/header").status_code == 404, "unknown run -> 404 (header)")
        check(client.get("/diff", params={"base": "caseA", "cand": "nope"}).status_code == 404,
              "/diff unknown cand -> 404")

    # -- missing-file path (503): app bound to a non-existent dataset -------
    client_missing = TestClient(create_app(os.path.join(tempfile.gettempdir(), "no_such_profile.h5")),
                                raise_server_exceptions=False)
    check(client_missing.get("/runs").status_code == 503, "missing dataset -> 503")

    print("\n" + "=" * 52)
    if problems == 0:
        print("TEST_API: ALL PASS (FastAPI backend serves the ProfilerResults contract)")
        return 0
    print(f"TEST_API: {problems} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
