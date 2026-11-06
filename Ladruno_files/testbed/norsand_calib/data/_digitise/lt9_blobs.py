"""Marker blobs of Lam & Tatsuoka Fig. 9 (list index = the marker ID used in run_lt9b.py).  Deterministic: same order each run."""
import pickle
from lt9_setup import *
P = lt9_panel()
bl = blobs_px(P, (0, 30), (-0.8, 8.0), rad=2)
pickle.dump(bl, open("lt9_blobs.pkl", "wb"))
print(len(bl), "blobs")
