"""Lab-curve schema (task item 3) and the synthetic writer.

CURVE CSV (one element test per file):
    columns  eps_a_pct, sr, eps_v_pct      (header row required; extra columns ignored)
      eps_a_pct  axial strain in %, compression positive (TE: as in the source, sign via eps_a_sign)
      sr         sigma1'/sigma3' (major/minor)
      eps_v_pct  volumetric strain in %; convention in the metadata (Tatsuoka's figures: dilation NEGATIVE)
    A row may leave sr or eps_v_pct empty (the two curves are digitised at different eps_a); each column is
    used where it is filled.
METADATA, either '# key: value' lines at the top of the CSV or a sidecar '<stem>.meta.json' (sidecar wins):
    required  id, kind (PS | TC | TE | PSU | TCU), sigma3_kPa (isotropic consolidation stress), e0 (void ratio
              at the start of shear; Tatsuoka quote e at 0.05 kgf/cm2 = 4.9 kPa, state which), source
              (paper + figure + page)
    optional  eps_v_convention ('dilation_negative' default | 'dilation_positive'),
              eps_a_sign (+1 default; -1 if the file's extension strain is positive),
              phi_peak_deg, eps_peak_pct (peak targets; default: read off the curve),
              weight (multiplies every residual of this test, default 1)
POINT-TEST CSV (phi_peak(e), eps_peak(e), e.g. Tatsuoka Figs. 9 and 22):
    columns  sigma3_kPa, e, phi_peak_deg, eps_peak_pct[, kind]   (eps_peak_pct may be empty; kind default PS)
"""
from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass, field

import numpy as np

REQUIRED_META = ("id", "kind", "sigma3_kPa", "e0", "source")


@dataclass
class LabCurve:
    meta: dict
    eps_a_sr: np.ndarray       # % compression positive, where sr is given
    sr: np.ndarray
    eps_a_ev: np.ndarray       # %, where eps_v is given
    eps_v: np.ndarray          # %, compression positive (dilation negative) after conversion
    path: str = ""

    @property
    def kind(self):
        return self.meta["kind"]

    @property
    def sigma3(self):
        return float(self.meta["sigma3_kPa"])

    @property
    def e0(self):
        return float(self.meta["e0"])

    def peak(self):
        """(sr_peak, eps_peak_pct): metadata if given, else the largest sr point of the curve."""
        if "phi_peak_deg" in self.meta:
            s = math.sin(math.radians(float(self.meta["phi_peak_deg"])))
            srp = (1.0 + s) / (1.0 - s)
        else:
            srp = float(np.max(self.sr))
        if "eps_peak_pct" in self.meta and self.meta["eps_peak_pct"] not in ("", None):
            ep = float(self.meta["eps_peak_pct"])
        else:
            ep = float(self.eps_a_sr[int(np.argmax(self.sr))])
        return srp, ep


@dataclass
class PointTest:
    sigma3: float
    e: float
    phi_peak_deg: float
    eps_peak_pct: float | None
    kind: str = "PS"
    source: str = ""


def _read_meta_lines(path):
    meta = {}
    with open(path, newline="") as f:
        for line in f:
            if not line.startswith("#"):
                break
            body = line[1:].strip()
            if ":" in body:
                k, v = body.split(":", 1)
                meta[k.strip()] = v.strip()
    return meta


def _coerce(meta):
    for k in ("sigma3_kPa", "e0", "phi_peak_deg", "eps_peak_pct", "weight", "eps_a_sign"):
        if k in meta and meta[k] not in ("", None):
            meta[k] = float(meta[k])
    return meta


def load_curve(path: str) -> LabCurve:
    meta = _read_meta_lines(path)
    side = os.path.splitext(path)[0] + ".meta.json"
    if os.path.exists(side):
        with open(side) as f:
            meta.update(json.load(f))
    meta = _coerce(meta)
    missing = [k for k in REQUIRED_META if k not in meta]
    if missing:
        raise ValueError(f"{path}: metadata missing {missing} (schema in harness/data.py)")
    if meta["kind"] not in ("PS", "TC", "TE", "PSU", "TCU"):
        raise ValueError(f"{path}: kind {meta['kind']!r} not in PS|TC|TE|PSU|TCU")
    rows = []
    with open(path, newline="") as f:
        rd = csv.DictReader(line for line in f if not line.startswith("#"))
        for r in rd:
            rows.append(r)

    def col(r, k):
        v = (r.get(k) or "").strip()
        return float(v) if v else float("nan")

    ea = np.array([col(r, "eps_a_pct") for r in rows]) * float(meta.get("eps_a_sign", 1.0))
    sr = np.array([col(r, "sr") for r in rows])
    ev = np.array([col(r, "eps_v_pct") for r in rows])
    conv = meta.get("eps_v_convention", "dilation_negative")
    if conv == "dilation_positive":
        ev = -ev
    elif conv != "dilation_negative":
        raise ValueError(f"{path}: eps_v_convention {conv!r}")
    ms, mv = ~np.isnan(sr) & ~np.isnan(ea), ~np.isnan(ev) & ~np.isnan(ea)
    o1, o2 = np.argsort(ea[ms], kind="stable"), np.argsort(ea[mv], kind="stable")
    return LabCurve(meta, ea[ms][o1], sr[ms][o1], ea[mv][o2], ev[mv][o2], path)


def write_curve(path: str, meta: dict, eps_a_pct, sr, eps_v_pct):
    """Writes the schema (meta as '# key: value' lines; eps_v in the dilation-negative convention)."""
    meta = dict(meta)
    meta.setdefault("eps_v_convention", "dilation_negative")
    with open(path, "w", newline="") as f:
        for k, v in meta.items():
            f.write(f"# {k}: {v}\n")
        w = csv.writer(f)
        w.writerow(["eps_a_pct", "sr", "eps_v_pct"])
        for a, s, e in zip(eps_a_pct, sr, eps_v_pct):
            w.writerow([repr(float(a)), repr(float(s)), repr(float(e))])


def load_points(path: str, source: str = "") -> list[PointTest]:
    out = []
    with open(path, newline="") as f:
        rd = csv.DictReader(line for line in f if not line.startswith("#"))
        for r in rd:
            ep = (r.get("eps_peak_pct") or "").strip()
            out.append(PointTest(float(r["sigma3_kPa"]), float(r["e"]), float(r["phi_peak_deg"]),
                                 float(ep) if ep else None, (r.get("kind") or "PS").strip() or "PS", source))
    return out
