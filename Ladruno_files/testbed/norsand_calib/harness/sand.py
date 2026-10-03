"""Constants of a sand that the P3 fit does NOT touch: the CSL, M, the elastic targets and p_a.

Plan 144 §3 P3 (re-aimed 2026-10-02): the CSL comes from DM04's published Toyoura Table 1 and is FIXED;
the elastic constants are fixed per the energy (energy.py maps the targets below onto the energy's own
constants); pi_i0 comes from the initial state (model.pi0_for). Units kPa, compression negative inside
the oracles; the targets here are written as positive magnitudes.

TOYOURA_PLACEHOLDER: Dafalias & Manzari (2004) J. Eng. Mech. 130(6):622-634, Table 1 (Toyoura), as recalled by
the harness author for the smoke gate (whose result does not depend on them). The data pack's transcription
(data/dm04/toyoura_table1.csv, read from the PDF journal p. 626, 2026-10-02) agrees on all seven values used here;
the placeholder's p_a = 101 kPa differs from the pack's p_at = 100 kPa (ASSUMED there: DM04 does not tabulate
p_at; data/README.md §4). Real fits use toyoura_dm04(), which reads the pack.
"""
from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass, field, asdict


@dataclass(frozen=True)
class ElasticTargets:
    """Physical elastic stiffness the energy plug must reproduce (DM04 form, DM04 eq. for G):
       G(p, e) = G0 p_a (2.97 - e)^2 / (1 + e) sqrt(p / p_a),   K = 2 (1 + nu) G / (3 (1 - 2 nu)).
    p is a positive magnitude (kPa)."""
    G0: float
    nu: float
    p_a: float
    n_exp: float = 0.5          # DM04: sqrt(p)

    def G(self, p: float, e: float) -> float:
        return self.G0 * self.p_a * (2.97 - e) ** 2 / (1.0 + e) * (p / self.p_a) ** self.n_exp

    def K(self, p: float, e: float) -> float:
        return 2.0 * (1.0 + self.nu) * self.G(p, e) / (3.0 * (1.0 - 2.0 * self.nu))


@dataclass(frozen=True)
class Sand:
    name: str
    M: float                     # critical stress ratio in compression (DM04 M_c -> NorSand M, sheet §15 'direct')
    c_ext: float                 # DM04 c = M_e/M_c (-> the natural pin for rho, sheet §15)
    e0: float                    # fork CSL e_c = e0 - lambda_c (p/p_a)^xi (sheet S.22)
    lambda_c: float
    xi: float
    p_a: float
    elastic: ElasticTargets
    source: str = ""
    notes: dict = field(default_factory=dict)

    def e_c(self, p: float) -> float:
        return self.e0 - self.lambda_c * (p / self.p_a) ** self.xi

    def as_dict(self) -> dict:
        d = asdict(self)
        return d


# Dafalias & Manzari (2004) Table 1, Toyoura: G0 125, nu 0.05, M 1.25, c 0.712, lambda_c 0.019, e0 0.934,
# xi 0.7, p_at 101 kPa (placeholder, see the module docstring).
TOYOURA_PLACEHOLDER = Sand(
    name="Toyoura (DM04 Table 1 as recalled, p_a 101: smoke placeholder; real fits use toyoura_dm04())",
    M=1.25, c_ext=0.712, e0=0.934, lambda_c=0.019, xi=0.7, p_a=101.0,
    elastic=ElasticTargets(G0=125.0, nu=0.05, p_a=101.0),
    source="Dafalias & Manzari 2004, JEM 130(6), Table 1 (recalled; matches data/dm04/toyoura_table1.csv except p_a)",
)


DATA_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))


def toyoura_dm04(p_a: float = 100.0, path: str | None = None) -> Sand:
    """Toyoura from the data pack's DM04 Table 1 transcription (data/dm04/toyoura_table1.csv; DM04 journal p. 626).
    p_a is NOT in DM04 Table 1: 100 kPa is the pack's stated assumption (data/README.md §4), passed explicitly."""
    path = path or os.path.join(DATA_DIR, "dm04", "toyoura_table1.csv")
    vals = {}
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            vals[r["constant"]] = float(r["value"])
    return Sand(name=f"Toyoura, DM04 Table 1 (data pack transcription), p_a {p_a:g} kPa (assumed)",
                M=vals["M"], c_ext=vals["c"], e0=vals["e0"], lambda_c=vals["lambda_c"], xi=vals["xi"], p_a=p_a,
                elastic=ElasticTargets(G0=vals["G0"], nu=vals["nu"], p_a=p_a),
                source="Dafalias & Manzari 2004 JEM 130(6) Table 1 p. 626, via data/dm04/toyoura_table1.csv")


def sand_from_dict(d: dict) -> Sand:
    """Build a Sand from a dict (e.g. the data agent's JSON). Keys as in the dataclass; 'elastic' a dict."""
    el = d["elastic"]
    return Sand(name=d["name"], M=float(d["M"]), c_ext=float(d["c_ext"]), e0=float(d["e0"]),
                lambda_c=float(d["lambda_c"]), xi=float(d["xi"]), p_a=float(d["p_a"]),
                elastic=ElasticTargets(G0=float(el["G0"]), nu=float(el["nu"]), p_a=float(el.get("p_a", d["p_a"])),
                                       n_exp=float(el.get("n_exp", 0.5))),
                source=d.get("source", ""), notes=d.get("notes", {}))


def phi_ps(sr: float) -> float:
    """Mobilised friction angle (deg) from sigma1'/sigma3' (Mohr-Coulomb, any loading mode)."""
    return math.degrees(math.asin((sr - 1.0) / (sr + 1.0)))
