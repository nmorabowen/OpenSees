"""TIMs ring states (WP-127 attachments) and the WP-128 smallest reproducer.

The ring CSVs live on the WP-127 branch at
Ladruno_implementation/_tims_2d_model_requests_2026-09-25/ring_points_b{8,16}.csv.
`load_ring_csv` reads the file from the checkout if present, else from git
(`git show <ref>:<path>`), else fails loudly.

SIGN: the sigma columns are INTERNAL, compression-positive (WP-127 finding A,
despite the attachment README); alpha, alpha_in, z are raw internal ratios.
"""
from __future__ import annotations

import csv
import io
import os
import subprocess

from .model import CAMPAIGN, Options, State

ATTACH = "Ladruno_implementation/_tims_2d_model_requests_2026-09-25"
GIT_REFS = ("origin/ladruno", "origin/wp/127-sanisand-replay-counters",
            "origin/wp/128-sanisand-ring-trace")


def _repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(here, "..", ".."))


def ring_csv_text(name):
    """name: 'b8' | 'b16' or a path."""
    if os.path.isfile(name):
        return open(name, newline="").read()
    rel = f"{ATTACH}/ring_points_{name}.csv"
    p = os.path.join(_repo_root(), rel)
    if os.path.isfile(p):
        return open(p, newline="").read()
    for ref in GIT_REFS:
        try:
            return subprocess.run(["git", "-C", _repo_root(), "show", f"{ref}:{rel}"],
                                  capture_output=True, text=True, check=True).stdout
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
    raise FileNotFoundError(f"ring CSV {name!r}: not in the checkout and not on {GIT_REFS}")


def load_ring_csv(name):
    rows = []
    for r in csv.DictReader(io.StringIO(ring_csv_text(name))):
        d = {k: float(v) for k, v in r.items()}
        d["element"], d["gp"] = int(d["element"]), int(d["gp"])
        for key in ("sigma", "alpha", "alpha_in", "z"):
            d[key] = [d[f"{key}_{i}"] for i in range(6)]
        d["mesh"] = name if name in ("b8", "b16") else os.path.basename(name)
        rows.append(d)
    return rows


def row_state(row):
    return State.from_voigt(row["sigma"], row["alpha"], row["z"], row["e"],
                            row["alpha_in"])


def probes(delta):
    """The documented WP-127 probes (compression positive, engineering shear,
    plane-strain admissible)."""
    return {"isoComp": [delta, delta, 0.0, 0.0, 0.0, 0.0],
            "shear": [0.0, 0.0, 0.0, delta, 0.0, 0.0]}


def reproducer_state(p_s=0.0101, e=None, P=CAMPAIGN):
    """WP-128's smallest reproducer: sigma = p_s I, alpha = alpha_in = z = 0."""
    import numpy as np
    z = np.zeros((3, 3))
    return State(p_s * np.eye(3), z.copy(), z.copy(), P.e_init if e is None else e, z.copy())


REPRODUCER_DEPS = [0.0, 1.0e-4, 0.0, 0.0, 0.0, 0.0]   # one plane-strain d eps_yy


def ring_variants(p_min=0.0101):
    """paper; uw_model (UW constitutive additions U1-U5, PAPER alpha_in rule,
    continuous moduli -- the oracle a corrected C++ integrator should hit);
    uw_rule (as uw_model + UW's once-per-increment alpha_in rule U6 and the h
    cap U7 -- shows mechanism G in the continuous model)."""
    uw_model = Options(d_factor=True, p_min=p_min, g_void_ratio="initial",
                       void_ratio_law="initial")
    return {"paper": Options(), "uw_model": uw_model,
            "uw_rule": uw_model.with_(alpha_in_rule="uw", h_cap=1.0e10)}
