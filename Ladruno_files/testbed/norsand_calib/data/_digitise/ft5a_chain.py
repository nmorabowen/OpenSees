"""Marker-chain run for Fukushima & Tatsuoka Fig. 5(a) sigma1'/sigma3' curves (reproduces outputs/ft5a_out.json['sr'] after the range / despike curation in ft5a_run.py).
(Originally run inline; reconstructed verbatim from the session -- the seeds are marker positions read from the zoomed scan at eps_a = 5 %.)"""
import pickle
from ft5a_setup import *
PF, P = filled_panel_ft5a()
bl = blobs_px(PF, (0.0, 15), (-0.5, 6.8))
seeds = {'0.1': (5, 6.10), '0.2': (5, 5.56), '0.5': (5, 5.31), '1.0': (5, 5.03), '2.0': (5, 4.64), '4.0': (5, 4.36)}
order = ['4.0', '2.0', '1.0', '0.1', '0.5', '0.2']
R = chain_family(PF, bl, seeds, (0, 15), (-0.5, 6.8), order=order)
pickle.dump(R, open("ft5a_chain.pkl", "wb"))
