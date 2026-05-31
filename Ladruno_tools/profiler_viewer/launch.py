#!/usr/bin/env python3
"""
One-click launcher for the Ladruno profiler viewer.

Stdlib-only bootstrapper (no third-party imports here): it provisions a private
virtualenv, installs the backend deps, builds the React frontend if a built
bundle is missing, then starts a single server that hosts BOTH the JSON API and
the UI and opens it in the browser. Cross-platform (Windows / Linux / macOS).

Usage:
    python launch.py [profile.h5] [--port N] [--no-browser] [--rebuild]

Normally invoked via the OS wrappers (Profiler_Viewer.bat / profiler_viewer.sh),
which just locate a Python and run this. If no profile.h5 is given it serves
whatever is in the current directory's profile.h5 (the UI shows a clear banner
when the file is missing, so it is fine to open it first and point it at a file
later by relaunching).
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
import venv
import webbrowser

HERE = os.path.dirname(os.path.abspath(__file__))
VENV_DIR = os.path.join(HERE, ".viewer_venv")
DIST_DIR = os.path.join(HERE, "frontend", "dist")
FRONTEND_DIR = os.path.join(HERE, "frontend")
REQS = os.path.join(HERE, "requirements.txt")
STAMP = os.path.join(VENV_DIR, ".deps_installed")  # marker so we don't reinstall every launch


def log(msg: str) -> None:
    print(f"[viewer] {msg}", flush=True)


def venv_python() -> str:
    """Path to the venv's python executable (platform-specific layout)."""
    if os.name == "nt":
        return os.path.join(VENV_DIR, "Scripts", "python.exe")
    return os.path.join(VENV_DIR, "bin", "python")


def ensure_venv() -> str:
    """Create the private venv if absent; return its python path."""
    py = venv_python()
    if not os.path.exists(py):
        log(f"creating virtualenv in {VENV_DIR} ...")
        venv.create(VENV_DIR, with_pip=True)
    return py


def ensure_backend_deps(py: str) -> None:
    """pip-install the backend requirements once (guarded by a stamp file)."""
    if os.path.exists(STAMP):
        return
    log("installing backend dependencies (fastapi, uvicorn, h5py, numpy) ...")
    subprocess.check_call([py, "-m", "pip", "install", "--quiet",
                           "--upgrade", "pip"])
    subprocess.check_call([py, "-m", "pip", "install", "--quiet", "-r", REQS])
    with open(STAMP, "w", encoding="utf-8") as f:
        f.write("ok\n")
    log("backend dependencies ready.")


def npm_cmd() -> str | None:
    """Locate npm (npm.cmd on Windows); None if Node isn't installed."""
    for name in (("npm.cmd", "npm") if os.name == "nt" else ("npm",)):
        path = shutil.which(name)
        if path:
            return path
    return None


def ensure_frontend(rebuild: bool) -> bool:
    """Build the React bundle into frontend/dist if missing. Returns True if a
    built UI is available, False if we must fall back to API-only."""
    index = os.path.join(DIST_DIR, "index.html")
    if os.path.exists(index) and not rebuild:
        return True

    npm = npm_cmd()
    if npm is None:
        log("WARNING: Node.js/npm not found — cannot build the UI.")
        log("         Install Node 18+ from https://nodejs.org and relaunch,")
        log("         or use the headless API at the printed URL (/docs).")
        return False

    if not os.path.exists(os.path.join(FRONTEND_DIR, "node_modules")):
        log("installing frontend packages (npm install, first run only) ...")
        subprocess.check_call([npm, "install"], cwd=FRONTEND_DIR)
    log("building the UI (npm run build) ...")
    subprocess.check_call([npm, "run", "build"], cwd=FRONTEND_DIR)
    return os.path.exists(index)


def main() -> int:
    ap = argparse.ArgumentParser(description="Launch the Ladruno profiler viewer.")
    ap.add_argument("h5", nargs="?", default=os.path.join(os.getcwd(), "profile.h5"),
                    help="path to the profiler profile.h5 (default: ./profile.h5)")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--rebuild", action="store_true",
                    help="force a fresh UI build even if frontend/dist exists")
    args = ap.parse_args()

    h5 = os.path.abspath(args.h5)
    if not os.path.exists(h5):
        log(f"NOTE: '{h5}' does not exist yet — the viewer will show a banner until")
        log("      a profile.h5 is created (run an analysis with profiler('report', ...)).")

    py = ensure_venv()
    ensure_backend_deps(py)
    have_ui = ensure_frontend(args.rebuild)

    url = f"http://{args.host}:{args.port}/"
    log(f"starting server at {url}  (Ctrl+C to stop)")
    if have_ui:
        log("the viewer UI is served at that URL.")
    else:
        log(f"API-only mode: open {url}docs for the raw endpoints.")

    # Launch the server as a child using the venv python so it sees the deps.
    # PROFILER_H5 tells profiler_api which dataset to serve; it mounts dist/ if present.
    env = dict(os.environ, PROFILER_H5=h5, LADRUNO_OPENSEES_QUIET="1")
    server = subprocess.Popen(
        [py, "-m", "uvicorn", "profiler_api:app",
         "--host", args.host, "--port", str(args.port), "--log-level", "warning"],
        cwd=HERE, env=env,
    )
    try:
        if not args.no_browser:
            time.sleep(1.5)          # give uvicorn a moment to bind the port
            webbrowser.open(url)
        server.wait()
    except KeyboardInterrupt:
        log("shutting down ...")
        server.terminate()
        try:
            server.wait(timeout=5)
        except subprocess.TimeoutExpired:
            server.kill()
    return server.returncode or 0


if __name__ == "__main__":
    sys.exit(main())
