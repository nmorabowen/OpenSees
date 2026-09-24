import os
import sys

# WP-105 / F12: this used to be a hardcoded absolute scratchpad path
# (C:\Users\nmb\...\release-build-634824e1f\dist\bin). Default now resolves to
# <fork root>/dist/bin (the build.bat REQUIRED layout, see CLAUDE.md); override
# with F12_DIST_BIN to point at a specific worktree's dist\bin.
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DIST = os.environ.get("F12_DIST_BIN", os.path.join(REPO_ROOT, "dist", "bin"))
sys.path.insert(0, DIST)
import opensees as ops
print("pyd:", ops.__file__)
print("build:", ops.ladrunoBuild())
