#!/usr/bin/env bash
# One-click launcher for the Ladruno profiler viewer (Linux / macOS).
# Run ./profiler_viewer.sh [profile.h5]. It finds Python, then hands off to
# launch.py (which provisions a venv, installs deps, builds the UI if needed,
# and opens the browser).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

PY=""
for c in python3 python; do
    if command -v "$c" >/dev/null 2>&1; then PY="$c"; break; fi
done
if [ -z "$PY" ]; then
    echo
    echo "  Python 3 was not found on PATH."
    echo "  Install Python 3.10+ (e.g. 'sudo apt install python3 python3-venv'"
    echo "  on Debian/Ubuntu, or from https://www.python.org), then run this again."
    echo
    exit 1
fi

exec "$PY" "$HERE/launch.py" "$@"
