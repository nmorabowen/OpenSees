"""Q2 plateau / slopes, Q3 B/16 localization fields, Q4 E_A vs E_B at TIMs' wall."""
import csv, json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

A = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(A, "..", "runs")
CK = os.path.join(A, "ck")
T = os.path.join(A, "tables")
B = 1.5
LEGS = ["E_A", "E_B", "E_D", "E_B16", "E_C2"]
col = dict(E_A="C3", E_B="C0", E_D="C2", E_B16="C1", E_C2="C4", DP="k")
lab = dict(E_A="E_A ModifiedEuler (IS 1)", E_B="E_B SAS-ME TolR 1e-4",
           E_D="E_D SAS-ME TolR 1e-3", E_B16="E_B16 SAS-ME, B/16 mesh",
           E_C2="E_C2 SAS-ME TanType 1 (running)", DP="ctrl DP 38 deg (local)")


def curve(leg):
    p = os.path.join(A, "ctrl_dp38_steps.csv") if leg == "DP" else os.path.join(R, leg, "steps.csv")
    S = list(csv.DictReader(open(p)))
    return (np.array([float(r["s_over_B"]) for r in S]), np.array([float(r["q_kPa"]) for r in S]),
            S)


def fit_slope(x, y, lo, hi):
    m = (x >= lo) & (x <= hi)
    if m.sum() < 3:
        return np.nan, int(m.sum())
    return np.polyfit(x[m], y[m], 1)[0], int(m.sum())


out = {}
C = {k: curve(k) for k in LEGS + ["DP"]}
# initial slope: fit over s/B in [5e-5, 5e-4] (step 1 carries the dead-load jump)
for k, (x, q, S) in C.items():
    k0, n0 = fit_slope(x, q, 5e-5, 5e-4)
    xe = x[-1]
    ends = {}
    for w in (0.001, 0.002, 0.005):
        ke, ne = fit_slope(x, q, xe - w, xe)
        ends[f"win{w}"] = dict(slope=ke, norm=ke / k0, n=ne)
    out[k] = dict(s_end=float(xe), q_end=float(q[-1]), q_max=float(q.max()), k0_kPa_per_sB=k0,
                  end_slopes=ends)
# slope history: sliding window of 0.002 s/B
fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
for k, (x, q, S) in C.items():
    ax[0].plot(x, q, "-", color=col[k], lw=1.3 if k != "DP" else 1.0, label=lab[k],
               ls="--" if k == "DP" else "-")
    ax[0].plot(x[-1], q[-1], "o", color=col[k], ms=4)
    k0 = out[k]["k0_kPa_per_sB"]
    xs, ss = [], []
    for xc in np.arange(0.002, x[-1] + 1e-9, 0.0005):
        kk, n = fit_slope(x, q, xc - 0.002, xc)
        if n >= 3:
            xs.append(xc); ss.append(kk / k0)
    ax[1].plot(xs, ss, color=col[k], label=lab[k], ls="--" if k == "DP" else "-")
    out[k]["slope_hist"] = [(round(a, 5), round(b, 4)) for a, b in zip(xs, ss)][::4]
ax[0].axvspan(0.026, 0.041, color="grey", alpha=.18, label="TIMs ME wall s/B 0.026-0.041")
ax[0].plot([0.116], [417.6], "s", color="m", label="TIMs PDMY limit 417.6 kPa @0.116")
ax[0].axhline(417.6, color="m", lw=.6, ls=":")
ax[0].set_xlabel("s/B"); ax[0].set_ylabel("q (kPa)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
ax[0].set_title("q-s, all legs (dots = last converged step, MODE FLOOR except C2)")
ax[1].set_yscale("log"); ax[1].set_xlabel("s/B (window end)"); ax[1].set_ylabel("dq/ds / initial slope")
ax[1].axvspan(0.026, 0.041, color="grey", alpha=.18)
ax[1].set_title("Tangent slope, 0.002 s/B sliding LSQ window, normalised by initial slope")
ax[1].legend(fontsize=7); ax[1].grid(alpha=.3, which="both")
fig.tight_layout(); fig.savefig(os.path.join(A, "qs_all_legs.png"), dpi=130)
ax[0].set_xlim(0, 0.06); ax[0].set_ylim(0, 1100); fig.savefig(os.path.join(A, "qs_all_legs_zoom.png"), dpi=130)


# ------------------------------------------------------------- Q3 fields ---
def load(leg, name):
    return np.load(os.path.join(CK, leg, "ckpt", name))


def gam(e):  # max in-plane engineering shear strain, from (exx, eyy, gxy)
    return np.sqrt((e[:, 0] - e[:, 1]) ** 2 + e[:, 2] ** 2)


def field_panel(ax, z, v, title, vmax, marks=None):
    m = (np.abs(z["gx"]) < 3.2) & (z["gy"] > -3.2)
    tri = mtri.Triangulation(z["gx"][m], z["gy"][m])
    cs = ax.tricontourf(tri, v[m], levels=np.linspace(0, vmax, 21), cmap="magma_r", extend="max")
    ax.plot([-0.75, 0.75], [0, 0], "c-", lw=4)
    if marks is not None:
        for (x, y, t) in marks:
            ax.plot(x, y, "x", color="lime", ms=9, mew=2)
    ax.set_aspect("equal"); ax.set_title(title, fontsize=8)
    return cs


def band_profiles(z, v, h, depths):
    """For each depth, the x-profile of v on the GP row nearest that depth (x>0 side),
    the peak location and the FWHM of the peak (in m and in element sizes h)."""
    res = []
    for y0 in depths:
        gy = z["gy"]
        yr = gy[np.argmin(np.abs(gy - y0))]
        m = (np.abs(gy - yr) < 1e-6) & (z["gx"] > 0) & (z["gx"] < 4)
        x = z["gx"][m]; y = v[m]; o = np.argsort(x); x, y = x[o], y[o]
        i = int(np.argmax(y)); pk = y[i]; base = np.median(y)
        half = base + 0.5 * (pk - base)
        lo = i
        while lo > 0 and y[lo] > half:
            lo -= 1
        hi = i
        while hi < len(y) - 1 and y[hi] > half:
            hi += 1
        # linear interpolation of the half crossings
        def cross(a, b):
            if y[a] == y[b]:
                return x[a]
            return x[a] + (half - y[a]) * (x[b] - x[a]) / (y[b] - y[a])
        xl = cross(lo, lo + 1) if y[lo] <= half else x[lo]
        xr = cross(hi - 1, hi) if y[hi] <= half else x[hi]
        res.append(dict(depth=float(yr), x_peak=float(x[i]), peak=float(pk), base=float(base),
                        peak_over_base=float(pk / base) if base > 0 else np.inf,
                        fwhm_m=float(xr - xl), fwhm_h=float((xr - xl) / h), xprof=x.tolist(), yprof=y.tolist()))
    return res


pairs = {  # interval increment + total push strain, matched s/B
    "E_B16 70->77 (s/B 0.01246->0.01352)": ("E_B16", "field_step00070.npz", "field_last_converged.npz", B / 16),
    "E_B 50->55 (s/B 0.01287->0.01404)": ("E_B", "field_step00050.npz", "field_step00055.npz", B / 8),
    "E_B 370->377 (s/B 0.05080->0.05084, at its wall)": ("E_B", "field_step00370.npz", "field_last_converged.npz", B / 8),
    "E_B 110->120 (s/B 0.0270->0.0294)": ("E_B", "field_step00110.npz", "field_step00120.npz", B / 8),
}
flo = list(csv.DictReader(open(os.path.join(T, "floor_refusers.csv"))))
fig, axs = plt.subplots(2, 4, figsize=(20, 9))
bands = {}
for j, (name, (leg, f0, f1, h)) in enumerate(pairs.items()):
    z0, z1, zf = load(leg, f0), load(leg, f1), load(leg, "field_flip.npz")
    dg = gam(z1["eps"] - z0["eps"])
    gt = gam(z1["eps"] - zf["eps"])
    marks = [(float(r["x"]), float(r["y"]), r["codes"]) for r in flo if r["leg"] == leg] \
        if "last" in f1 else None
    cs = field_panel(axs[0, j], z1, gt, f"{name}\nTOTAL push gamma_max (eps - eps_flip)",
                     np.percentile(gt, 99.7), marks)
    plt.colorbar(cs, ax=axs[0, j], shrink=.7)
    cs = field_panel(axs[1, j], z1, dg / max(float(z1["s_over_B"] - z0["s_over_B"]), 1e-12),
                     f"{name}\nINCREMENTAL d(gamma_max)/d(s/B) over the interval",
                     np.percentile(dg / max(float(z1["s_over_B"] - z0["s_over_B"]), 1e-12), 99.7), marks)
    plt.colorbar(cs, ax=axs[1, j], shrink=.7)
    depths = [-0.3, -0.6, -1.0, -1.5]
    bands[name] = dict(h=h, total=band_profiles(z1, gt, h, depths),
                       incr=band_profiles(z1, dg, h, depths),
                       max_gamma_total=float(gt.max()),
                       loc_max_total=(float(z1["gx"][np.argmax(gt)]), float(z1["gy"][np.argmax(gt)])))
fig.tight_layout(); fig.savefig(os.path.join(A, "fields_shear_strain.png"), dpi=110)

# profile plot: B/16 vs B/8 at matched s/B, several depths
fig, axs = plt.subplots(2, 4, figsize=(18, 7))
for i, kind in enumerate(("total", "incr")):
    for jd in range(4):
        a = axs[i, jd]
        for name, c in (("E_B16 70->77 (s/B 0.01246->0.01352)", "C1"), ("E_B 50->55 (s/B 0.01287->0.01404)", "C0")):
            b = bands[name][kind][jd]
            a.plot(b["xprof"], np.array(b["yprof"]) / b["peak"], ".-", color=c, ms=3,
                   label=f"{name.split(' ')[0]} FWHM {b['fwhm_m']:.2f} m = {b['fwhm_h']:.1f} h")
        a.set_title(f"{kind} gamma_max, GP row y = {b['depth']:.2f} m (normalised by peak)", fontsize=8)
        a.axvline(0.75, color="k", lw=.5); a.legend(fontsize=7); a.set_xlabel("x (m)")
fig.tight_layout(); fig.savefig(os.path.join(A, "band_profiles_B8_vs_B16.png"), dpi=110)
for v in bands.values():
    for kind in ("total", "incr"):
        for b in v[kind]:
            b.pop("xprof"); b.pop("yprof")
out["bands"] = bands

# E_B vs E_B16 curve at the same s/B
xb, qb, _ = C["E_B"]; x16, q16, _ = C["E_B16"]
cmp16 = []
for s in (0.002, 0.005, 0.008, 0.010, 0.012, 0.013, float(x16[-1])):
    cmp16.append(dict(s_over_B=s, q_B8=float(np.interp(s, xb, qb)), q_B16=float(np.interp(s, x16, q16)),
                      diff_pct=float(100 * (np.interp(s, x16, q16) / np.interp(s, xb, qb) - 1))))
out["B16_vs_B8_curve"] = cmp16


# refusals by region, per leg, over the whole run (log per-step lines are top-6 only,
# so use the census: cumulative refusals per point at the last converged step)
def region(x, y):
    ax_ = abs(x)
    if ax_ <= 0.75 and y > -0.2:
        return "under footing, top row"
    if ax_ <= 0.75:
        return "under footing (wedge/core)"
    if y > -0.2 and ax_ < 1.2:
        return "surface ring, edge (0.75<|x|<1.2)"
    if y > -0.2:
        return "surface ring, far (|x|>=1.2)"
    if ax_ < 1.2:
        return "edge fan below surface"
    return "other"


reg = {}
for leg in ("E_A", "E_B", "E_D", "E_B16"):
    z = load(leg, "field_last_converged.npz")
    st = z["stats"]
    col_ref = 9 if st.shape[1] == 17 else 15
    d = st[:, col_ref]
    rr = {}
    for k in np.where(d > 0)[0]:
        g = region(float(z["gx"][k]), float(z["gy"][k]))
        rr.setdefault(g, [0, 0]); rr[g][0] += int(d[k]); rr[g][1] += 1
    reg[leg] = dict(total=int(d.sum()), n_points=int((d > 0).sum()), by_region=rr,
                    n_points_total=int(len(d)))
    if col_ref == 15:
        names = {16: "startOutsideYield", 17: "startAlphaOutsideBounding", 19: "errorAtDTmin",
                 20: "loadingNonPosH", 21: "tensionAtDTmin", 22: "driftFailed", 23: "alphaOutsideAtDTmin",
                 24: "maxSubsteps"}
        reg[leg]["by_code"] = {v: int(st[:, i].sum()) for i, v in names.items() if st[:, i].sum() > 0}
out["refusals_by_region_whole_run"] = reg


# ------------------------------------------------ Q4 E_A vs E_B, s/B 0.026-0.030 ---
def win(leg, lo, hi):
    S = C[leg][2]
    W = [r for r in S if lo <= float(r["s_over_B"]) <= hi]
    return dict(steps=len(W), rungs={c: sum(1 for r in W if r["rung"] == str(i)) for i, c in enumerate("NLK")},
                iters_sum=sum(int(r["iters"]) for r in W),
                iters_median=float(np.median([int(r["iters"]) for r in W])) if W else None,
                failed_attempts=sum(int(r["fails_before"]) for r in W),
                refusals_or_caps=sum(int(r["cap_step"]) for r in W),
                forcedAtDTmin_or_refDTmin=sum(int(r["dtmin_step"]) for r in W),
                rejErr=sum(int(r["rejErr_step"]) for r in W),
                substeps=sum(int(r["sub_step_total"]) for r in W),
                ds_median=float(np.median([float(r["ds_m"]) for r in W])) if W else None,
                ds_min=min(float(r["ds_m"]) for r in W) if W else None,
                max_rho=max(float(r["max_rho_alpha"]) for r in W) if W else None,
                n_rho_gt_1_max=max(int(r["n_rho_gt_1"]) for r in W) if W else None,
                min_p=min(float(r["min_p"]) for r in W) if W else None,
                wall_h=sum(float(r["wall_step_s"]) for r in W) / 3600)


Q4 = {}
for lo, hi in ((0.020, 0.026), (0.026, 0.0293), (0.0293, 0.041)):
    Q4[f"{lo}-{hi}"] = {leg: win(leg, lo, hi) for leg in ("E_A", "E_B", "E_D")}


def fieldstats(leg, name):
    z = load(leg, name)
    rho, p = z["rho"], z["p"]
    return dict(s_over_B=float(z["s_over_B"]), step=int(z["step"]),
                n_rho_gt_1=int((rho > 1).sum()), n_rho_gt_0p99=int((rho > .99).sum()),
                n_rho_gt_1p01=int((rho > 1.01).sum()),
                rho_max=float(np.nanmax(rho)), rho_p99=float(np.nanpercentile(rho, 99)),
                p_min=float(p.min()), n_p_lt_5=int((p < 5).sum()), n_p_lt_10=int((p < 10).sum()),
                n_at_Pmin=int((p < 0.0102 * 1.0001).sum()),
                q_curve=float(np.interp(float(z["s_over_B"]), C[leg][0], C[leg][1])))


Q4["fields"] = {"E_A step130": fieldstats("E_A", "field_step00130.npz"),
                "E_B step105": fieldstats("E_B", "field_step00105.npz"),
                "E_B step110": fieldstats("E_B", "field_step00110.npz"),
                "E_A step185": fieldstats("E_A", "field_step00185.npz"),
                "E_A last(189)": fieldstats("E_A", "field_last_converged.npz"),
                "E_B step120": fieldstats("E_B", "field_step00120.npz"),
                "E_D step155": fieldstats("E_D", "field_step00155.npz")}
# q at matched s/B
Q4["q_matched"] = [dict(s_over_B=s, **{leg: float(np.interp(s, C[leg][0], C[leg][1])) for leg in ("E_A", "E_B", "E_D")})
                   for s in (0.010, 0.017, 0.020, 0.026, 0.028, 0.0292)]
# cumulative refusals/caps vs s/B for A and B through 0.030 (the per-point census at matched s/B)
zA, zB = load("E_A", "field_step00185.npz"), load("E_B", "field_step00120.npz")
Q4["cum_to_0.0292"] = dict(E_A_capHits=int(zA["stats"][:, 9].sum()), E_A_forcedAtDTmin=int(zA["stats"][:, 5].sum()),
                           E_A_capHit_points=int((zA["stats"][:, 9] > 0).sum()),
                           E_B_refusals=int(zB["stats"][:, 15].sum()), E_B_refDTmin=int(zB["stats"][:, 19].sum()),
                           E_B_refusal_points=int((zB["stats"][:, 15] > 0).sum()))
# where are A's forcedAtDTmin points (accepted with error above tolerance at dT_min)
fd = zA["stats"][:, 5]
top = np.argsort(-fd)[:10]
Q4["E_A_forcedAtDTmin_top10"] = [dict(x=round(float(zA["gx"][k]), 3), y=round(float(zA["gy"][k]), 3),
                                      n=int(fd[k]), p=round(float(zA["p"][k]), 2), rho=round(float(zA["rho"][k]), 4))
                                 for k in top]
out["Q4"] = Q4

# rho_alpha distributions at matched s/B, A vs B
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
for (leg, f, c) in (("E_A", "field_step00130.npz", "C3"), ("E_B", "field_step00105.npz", "C0"),
                    ("E_A", "field_step00185.npz", "darkred"), ("E_B", "field_step00120.npz", "navy")):
    z = load(leg, f)
    r = np.sort(z["rho"][np.isfinite(z["rho"])])[::-1]
    ax[0].semilogx(np.arange(1, len(r) + 1), r, color=c, label=f"{leg} step {int(z['step'])} s/B {float(z['s_over_B']):.4f}")
ax[0].axhline(1, color="k", lw=.6); ax[0].set_ylim(0.9, 1.04); ax[0].set_xlabel("rank"); ax[0].set_ylabel("rho_alpha")
ax[0].legend(fontsize=7); ax[0].set_title("rho_alpha, ranked (top of the distribution)")
for leg, c in (("E_A", "C3"), ("E_B", "C0"), ("E_D", "C2")):
    S = C[leg][2]
    x = [float(r["s_over_B"]) for r in S]
    ax[1].plot(x, [int(r["n_rho_gt_1"]) for r in S], color=c, label=f"{leg} n(rho>1)")
    ax[1].plot(x, np.cumsum([int(r["cap_step"]) for r in S]) / 10, color=c, ls=":", label=f"{leg} cum refusals/caps /10")
ax[1].axvspan(0.026, 0.041, color="grey", alpha=.18); ax[1].set_xlim(0, 0.052); ax[1].set_xlabel("s/B")
ax[1].legend(fontsize=7); ax[1].set_title("points with rho_alpha > 1, and cumulative refusals")
fig.tight_layout(); fig.savefig(os.path.join(A, "EA_vs_EB_rho.png"), dpi=130)

json.dump(out, open(os.path.join(T, "curves_fields_summary.json"), "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in out.items() if k not in ()}, indent=1, default=float)[:20000])
