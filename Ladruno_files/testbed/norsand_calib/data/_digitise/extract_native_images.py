"""Extract the embedded page bitmaps (1-bit, 2848 px wide, ~345 dpi) of the two scanned S&F papers into work/<tag>_native_<page>.png.

    python extract_native_images.py "<dir with the PDFs>"      (default: C:/Users/nmora/Dropbox/SOILS_rev/WP144_calibration)
tag lt = Lam_Tatsuoka_1988_SF28-1_89.pdf (journal page = PDF page + 88), tag ft = Fukushima_Tatsuoka_1984_SF24-4_30.pdf (journal page = PDF page + 29).
Page images are NOT redistributed in the repo (copyright); only the digitised numbers are.
"""
import os
import sys

import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1] if len(sys.argv) > 1 else "C:/Users/nmora/Dropbox/SOILS_rev/WP144_calibration"
OUT = os.environ.get("DIG_WORK") or os.path.join(HERE, "work")
os.makedirs(OUT, exist_ok=True)
for tag, name in (("lt", "Lam_Tatsuoka_1988_SF28-1_89.pdf"), ("ft", "Fukushima_Tatsuoka_1984_SF24-4_30.pdf")):
    doc = fitz.open(os.path.join(SRC, name))
    for i, pg in enumerate(doc):
        info = pg.get_image_info(xrefs=True)
        if not info:
            continue
        pix = fitz.Pixmap(doc, info[0]["xref"])
        if pix.n > 1:
            pix = fitz.Pixmap(fitz.csGRAY, pix)
        pix.save(os.path.join(OUT, f"{tag}_native_{i + 1:02d}.png"))
    print(tag, "done")
