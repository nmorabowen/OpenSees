"""The objective (task item 3): weighted residuals of the model against lab curves and point tests.

Per curve (kind, sigma3, e0 from its metadata; model = O2 driver, fixed axial increment):
  r_sr  = w_i (sr_model(eps_a_i) - sr_i) / SIGMA_SR           for every sr point with eps_a_i <= window end
  r_ev  = w_i (ev_model(eps_a_i) - ev_i) / SIGMA_EV           same for the eps_v points (%, compression +)
  r_phi = (phi_peak_model - phi_peak_data) / SIGMA_PHI         peak friction angle (deg)
  r_epk = (eps_peak_model - eps_peak_data) / SIGMA_EPK         axial strain at the peak (%)
  r_pen = REFUSAL_PENALTY * (1 - eps_reached / eps_needed)     0 when the driver completes (fixed length)
  window end = eps_peak_data + POST_PEAK_PCT ("up to and past the peak"); w_i = 1 up to the data peak,
  POST_PEAK_WEIGHT after it (post-peak lab response is localisation-affected; a single point is not).
Per point test (phi_peak(e), eps_peak(e), e.g. Tatsuoka Figs. 9/22): r_phi, r_epk (if given), r_pen.
Each curve's residuals are multiplied by its metadata 'weight' (default 1).
Model values at the data strains: linear interpolation on the model curve (|eps_a| axis). Model peak: vertex
of the parabola through the largest sr sample and its two neighbours (smooth in the parameters).

The default weights (Weights() below) are HARNESS CHOICES stated here, not data-derived:
  SIGMA_SR 0.10   ~ half the R_peak spread of two near-identical Tatsuoka specimens (PSL-24 e 0.714: 6.32 vs
                   Fig. 6a e 0.700: 6.52; data/tatsuoka1986/source/compare.md), >> digitisation 0.011
                   (source/tests_data.py docstring)
  SIGMA_EV 0.10 % the same order as the PSL-24 / Fig. 6a eps_v spread at eps_a 2-3 % (source/tests_data.py)
  SIGMA_PHI 0.5 deg, SIGMA_EPK 0.25 %: Fig. 9/22 read-off resolution and specimen scatter (owner to confirm)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict, field

import numpy as np

from .drivers import simulate
from .sand import phi_ps


@dataclass(frozen=True)
class Weights:
    sigma_sr: float = 0.10
    sigma_ev: float = 0.10
    sigma_phi_deg: float = 0.5
    sigma_eps_peak_pct: float = 0.25
    post_peak_pct: float = 2.0
    post_peak_weight: float = 0.5
    peak_terms: bool = True
    deps_pct: float = 0.05            # model axial increment, % (fixed: the residual is smooth in theta)
    margin_pct: float = 0.5           # model runs this far past the window end
    eps_min_pct: float = 0.0          # data points below this axial strain are left out (fit_tatsuoka uses 0.1:
    #                                   the (0, 1) first row is a TIMs prepend and the < 0.1 % branch carries
    #                                   plate seating, data/README.md §3.3)
    refusal_penalty: float = 100.0

    def as_dict(self):
        return asdict(self)


def model_peak(x: np.ndarray, sr: np.ndarray):
    """(sr_peak, x_peak, interior): parabola vertex through the max sample and its neighbours."""
    i = int(np.argmax(sr))
    if 0 < i < len(sr) - 1:
        x0, x1, x2 = x[i - 1:i + 2]
        y0, y1, y2 = sr[i - 1:i + 2]
        den = (x0 - x1) * (x0 - x2) * (x1 - x2)
        a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den
        b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / den
        c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1 + x0 * x1 * (x0 - x1) * y2) / den
        if a < 0.0:
            xv = -b / (2.0 * a)
            if x0 <= xv <= x2:
                return a * xv * xv + b * xv + c, xv, True
        return float(sr[i]), float(x[i]), True
    return float(sr[i]), float(x[i]), False


class Objective:
    def __init__(self, setup, curves=(), points=(), weights: Weights = Weights(), oracle: str = "O2"):
        self.setup = setup
        self.curves = list(curves)
        self.points = list(points)
        self.W = weights
        self.oracle = oracle
        self.n_evals = 0
        self.sim_seconds = 0.0

    # -- one curve -------------------------------------------------------------------------------
    def _curve_terms(self, theta, lc):
        W = self.W
        srp_d, epk_d = lc.peak()
        x_end = epk_d + W.post_peak_pct
        x_need = x_end + W.margin_pct
        n = int(math.ceil(x_need / W.deps_pct - 1e-9))
        # run to n * deps exactly: the model grid is k * deps for every test and every theta
        cv = simulate(self.setup, theta, lc.kind, lc.sigma3, lc.e0, n * W.deps_pct / 100.0, n, oracle=self.oracle)
        self.sim_seconds += cv.stats.get("seconds", 0.0)
        xm = np.abs(cv.eps_a_pct)
        srm, evm = cv.sr, cv.eps_v_pct
        wt = float(lc.meta.get("weight", 1.0))
        out = {}
        xs = np.abs(lc.eps_a_sr)
        ms = (xs <= x_end + 1e-12) & (xs >= W.eps_min_pct)
        ws = np.where(xs[ms] <= epk_d + 1e-12, 1.0, W.post_peak_weight)
        out["sr"] = wt * ws * (np.interp(xs[ms], xm, srm) - lc.sr[ms]) / W.sigma_sr
        xv = np.abs(lc.eps_a_ev)
        mv = (xv <= x_end + 1e-12) & (xv >= W.eps_min_pct)
        wv = np.where(xv[mv] <= epk_d + 1e-12, 1.0, W.post_peak_weight)
        out["ev"] = wt * wv * (np.interp(xv[mv], xm, evm) - lc.eps_v[mv]) / W.sigma_ev
        srp_m, xpk_m, interior = model_peak(xm, srm)
        if W.peak_terms:
            out["peak"] = wt * np.array([(phi_ps(srp_m) - phi_ps(srp_d)) / W.sigma_phi_deg,
                                         (xpk_m - epk_d) / W.sigma_eps_peak_pct])
        frac = min(1.0, xm[-1] / x_need) if len(xm) > 1 else 0.0
        out["pen"] = np.array([W.refusal_penalty * (1.0 - frac) if not cv.complete else 0.0])
        info = dict(id=lc.meta["id"], status=cv.status, phi_model=phi_ps(srp_m), phi_data=phi_ps(srp_d),
                    eps_peak_model=xpk_m, eps_peak_data=epk_d, peak_interior=interior,
                    rms_sr=float(np.sqrt(np.mean(out["sr"] ** 2))) if out["sr"].size else 0.0,
                    rms_ev=float(np.sqrt(np.mean(out["ev"] ** 2))) if out["ev"].size else 0.0,
                    seconds=cv.stats.get("seconds", 0.0))
        return out, info, cv

    def _point_terms(self, theta, pt):
        W = self.W
        ep = pt.eps_peak_pct if pt.eps_peak_pct is not None else 3.0
        x_need = 2.0 * ep + 1.0
        n = int(math.ceil(x_need / W.deps_pct - 1e-9))
        cv = simulate(self.setup, theta, pt.kind, pt.sigma3, pt.e, n * W.deps_pct / 100.0, n, oracle=self.oracle)
        self.sim_seconds += cv.stats.get("seconds", 0.0)
        xm = np.abs(cv.eps_a_pct)
        srp_m, xpk_m, interior = model_peak(xm, cv.sr)
        r = [(phi_ps(srp_m) - pt.phi_peak_deg) / W.sigma_phi_deg]
        if pt.eps_peak_pct is not None:
            r.append((xpk_m - pt.eps_peak_pct) / W.sigma_eps_peak_pct)
        frac = min(1.0, xm[-1] / x_need) if len(xm) > 1 else 0.0
        r.append(W.refusal_penalty * (1.0 - frac) if not cv.complete else 0.0)
        info = dict(sigma3=pt.sigma3, e=pt.e, status=cv.status, phi_model=phi_ps(srp_m), phi_data=pt.phi_peak_deg,
                    eps_peak_model=xpk_m, eps_peak_data=pt.eps_peak_pct, peak_interior=interior)
        return np.array(r), info

    def noise_variance(self) -> np.ndarray:
        """Variance of each residual (same order as residuals()) when the DATA carry independent Gaussian noise of
        standard deviation sigma_* (the Weights): r = w (model - data)/sigma, so var(r_i) = w_i^2, the curve weight
        included; the peak terms have unit weight, the refusal-penalty entries none (0). Uses the same masks as
        _curve_terms (windows come from the data peak, so a noisy peak target moves the window slightly: the
        variance is evaluated for the objective's own curves)."""
        W = self.W
        parts = []
        for lc in self.curves:
            srp_d, epk_d = lc.peak()
            x_end = epk_d + W.post_peak_pct
            wt = float(lc.meta.get("weight", 1.0))
            xs = np.abs(lc.eps_a_sr)
            ms = (xs <= x_end + 1e-12) & (xs >= W.eps_min_pct)
            ws = np.where(xs[ms] <= epk_d + 1e-12, 1.0, W.post_peak_weight)
            xv = np.abs(lc.eps_a_ev)
            mv = (xv <= x_end + 1e-12) & (xv >= W.eps_min_pct)
            wv = np.where(xv[mv] <= epk_d + 1e-12, 1.0, W.post_peak_weight)
            parts += [(wt * ws) ** 2, (wt * wv) ** 2]
            if W.peak_terms:
                parts.append(np.full(2, wt ** 2))
            parts.append(np.zeros(1))
        for pt in self.points:
            parts.append(np.array([1.0] + ([1.0] if pt.eps_peak_pct is not None else []) + [0.0]))
        return np.concatenate(parts) if parts else np.zeros(0)

    # -- public -----------------------------------------------------------------------------------
    def residuals(self, theta: dict) -> np.ndarray:
        self.n_evals += 1
        parts = []
        for lc in self.curves:
            out, _, _ = self._curve_terms(theta, lc)
            parts += [out["sr"], out["ev"]] + ([out["peak"]] if "peak" in out else []) + [out["pen"]]
        for pt in self.points:
            r, _ = self._point_terms(theta, pt)
            parts.append(r)
        return np.concatenate(parts) if parts else np.zeros(0)

    def breakdown(self, theta: dict) -> dict:
        curves, points = [], []
        total = 0.0
        for lc in self.curves:
            out, info, _ = self._curve_terms(theta, lc)
            info["cost_part"] = 0.5 * float(sum(np.sum(v ** 2) for v in out.values()))
            total += info["cost_part"]
            curves.append(info)
        for pt in self.points:
            r, info = self._point_terms(theta, pt)
            info["cost_part"] = 0.5 * float(np.sum(r ** 2))
            total += info["cost_part"]
            points.append(info)
        return dict(cost=total, curves=curves, points=points, weights=self.W.as_dict())
