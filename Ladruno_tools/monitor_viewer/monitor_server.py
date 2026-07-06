"""
monitor_server.py  —  live web dashboard for a Ladruno Analysis-monitor sink.

A thin FastAPI layer over MonitorReader that serves a self-contained HTML page
(monitor_page.html) which polls for new frames and draws them as live line
charts — the "watch it in a browser while the run streams" experience.

    python monitor_server.py sink.h5 [--host 127.0.0.1] [--port 8030]
    # or:  MONITOR_H5=sink.h5 uvicorn monitor_server:app

Endpoints:
    GET /                     the dashboard page
    GET /health              {status, file, nframes}
    GET /api/meta            columns + counts + units
    GET /api/frames?since=N  incremental frames [N:n)  (the poll endpoint)

Reads the sink per request (SWMR refresh), so it tails a file another process
is still writing and never holds a long-lived lock.
"""

from __future__ import annotations

import argparse
import os

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse

from monitor_reader import MonitorReader

_HERE = os.path.dirname(os.path.abspath(__file__))
_PAGE = os.path.join(_HERE, "monitor_page.html")


def create_app(sink_path: str) -> FastAPI:
    app = FastAPI(title="Ladruno Monitor Viewer")
    app.add_middleware(
        CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"],
    )
    reader = MonitorReader(sink_path)

    @app.get("/health")
    def health():
        try:
            m = reader.meta()
            return {"status": "ok", "file": m["path"], "nframes": m["nframes"]}
        except FileNotFoundError:
            return JSONResponse(status_code=503,
                                content={"status": "waiting", "file": sink_path,
                                         "detail": "sink not created yet"})
        except Exception as e:
            return JSONResponse(status_code=503,
                                content={"status": "error", "detail": str(e)})

    @app.get("/api/meta")
    def meta():
        try:
            return reader.meta()
        except FileNotFoundError:
            raise HTTPException(503, f"sink not created yet: {sink_path}")
        except ValueError as e:
            raise HTTPException(400, str(e))

    @app.get("/api/frames")
    def frames(since: int = Query(0, ge=0)):
        try:
            return reader.frames_since(since)
        except FileNotFoundError:
            raise HTTPException(503, f"sink not created yet: {sink_path}")
        except Exception as e:
            raise HTTPException(500, str(e))

    @app.get("/")
    def index():
        return FileResponse(_PAGE)

    return app


# Module-level app for `uvicorn monitor_server:app` (sink via $MONITOR_H5).
_env = os.environ.get("MONITOR_H5")
app = create_app(_env) if _env else None


def main(argv=None):
    ap = argparse.ArgumentParser(description="Serve the live monitor dashboard.")
    ap.add_argument("sink", help="path to the monitor SWMR-HDF5 sink")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8030)
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args(argv)

    import uvicorn
    application = create_app(args.sink)
    url = f"http://{args.host}:{args.port}"
    print(f"[monitor-viewer] serving {os.path.abspath(args.sink)} at {url}")
    if not args.no_browser:
        try:
            import webbrowser
            webbrowser.open(url)
        except Exception:
            pass
    uvicorn.run(application, host=args.host, port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
