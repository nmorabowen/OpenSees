"""WP-150 T5 control: the fork's WP-133 PDMY03 STAND-IN (NOT TIMs' PDMY01 33 deg calibration, which the fork does
not hold) in drained plane-strain compression at p0 = 10, 50, 150, 500 kPa, on one `quad` (tests/wp133_pdmy03_deck.py).
Reports phi'_peak, the axial strain at peak and the peak dilation angle, as t5_element_physics.py does for SANISAND.

Run with the engine under test (python -S, manual paths: the boot .pth preloads a stale pyd otherwise):
    <py3.12> -S t5_pdmy_control.py <dist/bin> [eps_max]
"""
import json, math, os, sys

BIN = os.path.abspath(sys.argv[1])
EPS_MAX = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
if sys.platform == "win32":
    os.add_dll_directory(BIN)
    sys.path.append(r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.join(ROOT, "tests"))
import numpy as np  # noqa: E402
import opensees as ops  # noqa: E402
assert os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) == os.path.normcase(BIN), ops.__file__
import wp133_pdmy03_deck as deck  # noqa: E402

rows, paths = [], {}
for p0 in (10.0, 50.0, 150.0, 500.0):
    deck.P0 = p0
    nstep = int(round(EPS_MAX / abs(deck.DY)))
    res = deck.run_quad(ops, 1, nstep=nstep, materials=lambda o: deck.build_material(o, 1, tail=deck.TAIL))
    S = np.array([s[:3] for s in res["stress"]])           # GP 1: sxx syy sxy (tension positive)
    E = np.array([e[:3] for e in res["strain"]])
    s1, s3 = -S[:, 1], -S[:, 0]                              # vertical = major, lateral = minor (compression +)
    phi = np.degrees(np.arcsin(np.clip((s1 - s3) / (s1 + s3), -1, 1)))
    ea, el = -E[:, 1], -E[:, 0]
    ev = ea + el
    de1, dev, de3 = np.gradient(ea), np.gradient(ev), np.gradient(el)
    with np.errstate(invalid="ignore", divide="ignore"):
        dil = np.degrees(np.arcsin(np.clip(-dev / (de1 - de3), -1, 1)))
    ip = int(np.argmax(phi))
    r = dict(p0=p0, n=int(res["n"]), phi_peak=float(phi[ip]), eps_a_peak=float(ea[ip]), phi_end=float(phi[-1]),
             eps_a_end=float(ea[-1]), dil_max=float(np.nanmax(dil)), ev_end=float(ev[-1]))
    rows.append(r)
    paths[int(p0)] = dict(eps_a=ea.tolist(), ev=ev.tolist(), phi=phi.tolist())
    print(f"PDMY03 stand-in PS p0={p0:5.0f}: steps {r['n']}  phi_peak {r['phi_peak']:.1f} @ eps_a {r['eps_a_peak']:.4f}  "
          f"phi_end {r['phi_end']:.1f} @ {r['eps_a_end']:.3f}  psi_max {r['dil_max']:.1f} deg  ev_end {r['ev_end']:+.4f}",
          flush=True)

L = ["| PDMY03 stand-in, PS | p0 kPa | φ′_peak ° | ε_a at peak | φ′ end ° (ε_a) | ψ_max ° | Bolton 0.8·ψ_max vs Δφ(peak − end) |",
     "|---|---|---|---|---|---|---|"]
for r in rows:
    L.append(f"| | {r['p0']:.0f} | {r['phi_peak']:.1f} | {r['eps_a_peak']:.4f} | {r['phi_end']:.1f} ({r['eps_a_end']:.3f}) | "
             f"{r['dil_max']:.1f} | {0.8 * r['dil_max']:.1f} vs {r['phi_peak'] - r['phi_end']:.1f} |")
txt = ("WP-133 PDMY03 stand-in (fork's, NOT TIMs' PDMY01 33°): φ 40°, PT 26°, G 1.3e5, B 2.6e5 at refP 101, d 0.5.\n\n"
       + "\n".join(L) + "\n")
print(txt)
open(os.path.join(HERE, "out_t5_pdmy03_standin.md"), "w", encoding="utf-8").write(txt)
json.dump(dict(rows=rows, paths=paths), open(os.path.join(HERE, "out_t5_pdmy03_standin.json"), "w"))
