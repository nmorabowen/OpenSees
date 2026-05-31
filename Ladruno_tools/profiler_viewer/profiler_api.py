"""
profiler_api.py  —  FastAPI backend over ProfilerResults (profiler phase P7).

Thin HTTP layer that exposes the headless ProfilerResults loader as JSON
endpoints for the React frontend (P8) and for ad-hoc curl/Jupyter use. It owns no
analysis logic of its own — every endpoint just forwards to a ProfilerResults
method, so the wire format is exactly what profiler_results.py already documents.

Run against one profile.h5 (one immutable group per run inside it):

    # module-level app reads $PROFILER_H5 (default ./profile.h5)
    PROFILER_H5=profile_sample.h5 uvicorn profiler_api:app --reload

    # or the CLI entry point with explicit flags
    python profiler_api.py --file profile_sample.h5 --port 8000

Then GET http://127.0.0.1:8000/docs for the live Swagger UI, or:

    GET /health                      liveness + file + run list
    GET /runs                        manifest rows           (picker)
    GET /runs/{run}/header           full run header
    GET /runs/{run}/rollup           normalized rollup tree  (flame graph)
    GET /runs/{run}/series           per-step arrays         (time history; 404 if none)
    GET /runs/{run}/memory           end-of-run memory snapshot + census
    GET /diff?base=&cand=            per-node wall-time delta (prove a fix)

The HDF5 file is opened per request (open -> read -> close) rather than held
open, so a run appended by a later `profiler('report', ...)` shows up without a
server restart, and there is no long-lived file lock (HDF5_USE_FILE_LOCKING).
"""

from __future__ import annotations

import os
from contextlib import contextmanager

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware

from profiler_results import ProfilerResults

DEFAULT_H5 = "profile.h5"


def create_app(h5_path: str) -> FastAPI:
    """Build a FastAPI app bound to a single profiler HDF5 dataset."""
    app = FastAPI(
        title="OpenSees/Ladruno profiler backend",
        version="1.0",
        summary="HTTP view over a profiler profile.h5 (P7).",
    )

    # Local dev/analysis tool: the React frontend runs on a different localhost
    # port, so allow any origin. This server should never be exposed publicly.
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_methods=["GET"],
        allow_headers=["*"],
    )

    @contextmanager
    def _open():
        """Open the dataset for one request; 503 if the file is missing/invalid."""
        if not os.path.exists(h5_path):
            raise HTTPException(status_code=503,
                                detail=f"profiler dataset not found: {h5_path}")
        try:
            pr = ProfilerResults(h5_path)
        except (OSError, ValueError) as e:
            raise HTTPException(status_code=503,
                                detail=f"cannot open {h5_path}: {e}")
        try:
            yield pr
        finally:
            pr.close()

    def _require_run(pr: ProfilerResults, run: str) -> None:
        if run not in pr.runs():
            raise HTTPException(status_code=404, detail=f"no such run: {run!r}")

    @app.get("/health")
    def health() -> dict:
        with _open() as pr:
            return {"status": "ok", "file": os.path.abspath(h5_path),
                    "schema_version": pr.schema_version, "runs": pr.runs()}

    @app.get("/runs")
    def runs() -> list[dict]:
        with _open() as pr:
            return pr.manifest()

    @app.get("/runs/{run}/header")
    def header(run: str) -> dict:
        with _open() as pr:
            _require_run(pr, run)
            return pr.header(run)

    @app.get("/runs/{run}/rollup")
    def rollup(run: str) -> dict:
        with _open() as pr:
            _require_run(pr, run)
            return pr.rollup(run)

    @app.get("/runs/{run}/series")
    def series(run: str) -> dict:
        with _open() as pr:
            _require_run(pr, run)
            s = pr.series(run)
            if s is None:
                raise HTTPException(
                    status_code=404,
                    detail=f"run {run!r} has no per-step series (run without -perStep)")
            return s

    @app.get("/runs/{run}/memory")
    def memory(run: str) -> dict:
        with _open() as pr:
            _require_run(pr, run)
            return pr.memory(run)

    @app.get("/diff")
    def diff(base: str = Query(...), cand: str = Query(...)) -> dict:
        with _open() as pr:
            _require_run(pr, base)
            _require_run(pr, cand)
            return {"base": base, "cand": cand, "rows": pr.diff(base, cand)}

    # Serve the built React frontend (frontend/dist) at the site root, when it
    # exists, so a single process hosts both the API and the UI (the one-click
    # launcher builds dist/ then starts this). html=True makes "/" return
    # index.html; SPA routing has no server routes so the API paths above win.
    # When dist/ is absent (dev: use the Vite server on :5173) this is skipped.
    dist = os.path.join(os.path.dirname(os.path.abspath(__file__)), "frontend", "dist")
    if os.path.isdir(dist):
        from fastapi.staticfiles import StaticFiles
        app.mount("/", StaticFiles(directory=dist, html=True), name="frontend")

    return app


# Module-level app for `uvicorn profiler_api:app`. The dataset path comes from
# $PROFILER_H5 so the server can be pointed at any profile.h5 without code edits.
app = create_app(os.environ.get("PROFILER_H5", DEFAULT_H5))


def main() -> None:
    import argparse

    import uvicorn

    ap = argparse.ArgumentParser(description="Serve a profiler profile.h5 over HTTP (P7).")
    ap.add_argument("--file", default=os.environ.get("PROFILER_H5", DEFAULT_H5),
                    help="path to the profiler HDF5 dataset (default: $PROFILER_H5 or profile.h5)")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8000)
    args = ap.parse_args()

    uvicorn.run(create_app(args.file), host=args.host, port=args.port)


if __name__ == "__main__":
    main()
