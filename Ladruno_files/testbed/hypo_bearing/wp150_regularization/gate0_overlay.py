"""WP-150 GATE 0 (a): overlay the exact oracle (and/or the C++) on Dafalias & Manzari (2004) Figs. 5-9, using the
reader's OWN licensed copy of the paper. The rendered page images and overlays are NOT committed (ASCE copyright).
    python gate0_overlay.py <dafalias2004.pdf> <out_gate0_oracle.json | merged json> <set[,set]> <outdir>
Axis calibration (CAL) was measured on 200-dpi renders from equally spaced dotted gridlines (median-band darkness
profiles), <= 1 px, and each axis is linear between two gridlines of known value.
"""
import json, sys
import numpy as np
from PIL import Image, ImageDraw
PANELS = {
    "F5a": (8, (100, 70, 830, 470)), "F5b": (8, (100, 505, 830, 1000)),
    "F6a": (8, (880, 70, 1700, 520)), "F6b": (8, (880, 565, 1700, 955)),
    "F7a": (8, (880, 1025, 1700, 1480)),
    "F8a": (9, (100, 60, 820, 490)), "F8b": (9, (100, 540, 820, 960)),
    "F9a": (9, (930, 60, 1700, 490)), "F9b": (9, (930, 540, 1700, 960)),
}
# (x: pix0, value0, pix1, value1), (y: pix0, value0, pix1, value1), x-quantity, tests
CAL = {
    "F5a": ((131, 0, 680, 3000), (370, 0, 37, 5000), "p", "F5"),
    "F5b": ((132, 0, 690, 30), (447, 0, 109.5, 5000), "ea", "F5"),
    "F6a": ((111, 0, 665, 3000), (425, 0, 89, 2500), "p", "F6"),
    "F6b": ((125, 0, 676, 30), (362, 0, 28, 2500), "ea", "F6"),
    "F7a": ((114, 0, 666, 2000), (428, 0, 93, 1200), "p", "F7"),
    "F8a": ((201, 0.80, 628.5, 0.95), (389, 0, 45, 1600), "e", "F8"),
    "F8b": ((132, 0, 690, 30), (344.3, 0, 6, 1600), "ea", "F8"),
    "F9a": ((63.5, 0.80, 622, 1.00), (383, 0, 45, 350), "e", "F9"),
    "F9b": ((58, 0, 626, 30), (347, 0, 2, 350), "ea", "F9"),
}
import os, fitz
pdf, srcp, setsarg, outdir = sys.argv[1:5]
os.makedirs(outdir, exist_ok=True)
doc = fitz.open(pdf)
for pno in (8, 9):
    doc[pno - 1].get_pixmap(dpi=200).save(os.path.join(outdir, f"page{pno}.png"))
src = json.load(open(srcp))
sets = setsarg.split(",")                 # e.g. paper  or  paper,cxx_off
colors = {"paper": (220, 0, 0), "uw_model": (0, 160, 0), "cxx_off": (0, 90, 255), "cxx_on": (255, 140, 0)}
for name, (page, box) in PANELS.items():
    (x0p, x0v, x1p, x1v), (y0p, y0v, y1p, y1v), xq, fig = CAL[name]
    im = Image.open(os.path.join(outdir, f"page{page}.png")).convert("RGB").crop(box)
    d = ImageDraw.Draw(im)
    def X(v): return x0p + (v - x0v) * (x1p - x0p) / (x1v - x0v)
    def Y(v): return y0p + (v - y0v) * (y1p - y0p) / (y1v - y0v)
    for s in sets:
        for label, rec in src[s].items():
            if not label.startswith(fig + "_"):
                continue
            xs = {"p": rec["p"], "e": rec["e"], "ea": [100 * v for v in rec["eps_a"]]}[xq]
            pts = [(X(a), Y(b)) for a, b in zip(xs, rec["q"])]
            d.line(pts, fill=colors[s], width=2)
    im.save(os.path.join(outdir, f"overlay_{name}_{'_'.join(sets)}.png"))
print("ok")
