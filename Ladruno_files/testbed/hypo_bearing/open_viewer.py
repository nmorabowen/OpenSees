"""Open the apeGmsh DESKTOP viewer on the bearing benchmark results.

Reads `bearing_view.ladruno` (written by viewer_run.py). The `.ladruno` route is
self-sufficient — the recorder file carries its own MODEL/NODES/ELEMENTS — so no
`model.h5` and no deck-to-FEMData tag reconciliation is needed.

Blocking VTK+Qt window: run it from a terminal (or detached), NOT from a
notebook, where it would take the kernel down. In a notebook use
`results.show_web()` or `results.viewer(blocking=False)` instead.

Run with the apeGmsh env, which has the viewer extra:
    C:/Users/nmb/venv/opensees_env/Scripts/python.exe open_viewer.py
"""
import os
import sys

from apeGmsh import Results

HERE = os.path.dirname(os.path.abspath(__file__))
FILE = os.path.join(HERE, sys.argv[1] if len(sys.argv) > 1
                    else "bearing_view.ladruno")

if not os.path.exists(FILE):
    raise SystemExit(f"{FILE} missing — run viewer_run.py first")

print(f"loading {os.path.basename(FILE)} "
      f"({os.path.getsize(FILE) / 1e6:.1f} MB) ...", flush=True)
results = Results.from_ladruno(FILE)
print("model:", results.model.fem.info, flush=True)
print("\nBenchmark: B = 2 m square smooth footing, 20 x 20 x 8 m graded domain,",
      flush=True)
print("2816 LadrunoUP H8, saturated PDMY (phi = 33 deg nominal), -geom hypo,",
      flush=True)
print("drained displacement-controlled push to s/B = 3 % (q = 1337 kPa).",
      flush=True)
print("\nFields recorded: nodal displacement; per-Gauss-point EFFECTIVE stress",
      flush=True)
print("and porePressure. 24 frames.", flush=True)
print("\nNote for reading the stress fields: the material's failure surface is a",
      flush=True)
print("Drucker-Prager CONE and it mobilises sin(phi) = 0.770 at failure, not",
      flush=True)
print("sin(33 deg) = 0.545 -- see vault note 17. Peak mobilisation anywhere in",
      flush=True)
print("this run is ~0.69, i.e. ~89 % of strength, and NO element is within 5 %",
      flush=True)
print("of failure.\n", flush=True)
print("opening the desktop viewer (close the window to exit) ...", flush=True)

results.viewer()
