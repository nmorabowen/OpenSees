"""WP-105 / F12 phase (b)/(c): the ADR-95/CP1 bearing deck with IntScheme 2.

`sanisand_tau0_band.py` hard-codes `INT_SCHEME = 1` as a module constant with no
CLI or env seam, and it is read at CALL time (`:696`, `:1090`). This wrapper
loads that driver as a module, sets the constant, and calls its own `main()` --
so NO file in the repo is edited and the driver's own provenance checks, controls
and JSON record still run.
"""
from __future__ import annotations

import importlib.util
import os
import sys

# WP-105 / F12: was a hardcoded absolute scratchpad path. Default now
# resolves to the sibling driver two levels up from this script
# (Ladruno_files/testbed/hypo_bearing/adr92_f12/ -> .../hypo_bearing/
# sanisand_tau0_band.py); override with F12_DRIVER.
_HYPO_BEARING = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DRIVER = os.environ.get("F12_DRIVER", os.path.join(_HYPO_BEARING, "sanisand_tau0_band.py"))
# Default resolves to <fork root>/dist/bin (the build.bat REQUIRED
# layout, see CLAUDE.md); override with F12_DIST_BIN.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DIST = os.environ.get("F12_DIST_BIN", os.path.join(_REPO_ROOT, "dist", "bin"))

os.environ["LADRUNO_DIST_BIN"] = DIST
os.environ["LADRUNO_A2_EXPECT_BUILD"] = "634824e1fbcf802bf27c7fbd29c113e5a2d63cb6"
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

scheme = int(sys.argv[1])
argv = sys.argv[2:]

spec = importlib.util.spec_from_file_location("tau0band", DRIVER)
mod = importlib.util.module_from_spec(spec)
sys.modules["tau0band"] = mod
spec.loader.exec_module(mod)

mod.INT_SCHEME = scheme
print("@@F12 INT_SCHEME forced to %d (driver default was 1)" % scheme, flush=True)
raise SystemExit(mod.main(argv))
