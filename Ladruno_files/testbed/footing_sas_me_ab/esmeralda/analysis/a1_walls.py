"""Q1 + Q5: why each leg hit the floor. Last-20-step tables, final-ladder census
diff (codes + locations), step-size history, Krylov-10x fractions."""
import csv, json, os, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(A, "..", "runs")
CK = os.path.join(A, "ck")
sys.path.insert(0, os.path.join(A, ".."))
import rung_fail  # noqa

LEGS = ["E_A", "E_B", "E_D", "E_B16", "E_C2"]
SAS = ["updates", "elastic", "substeps", "accepted", "rejectedErr", "rejectedLowP",
       "rejectedNonPosH", "rejectedDrift", "rejectedAlpha", "elasticStages",
       "driftCorrections", "alphaInReseats", "hBrackets", "alphaProjected",
       "intersectFail", "refusals", "refStartF", "refStartAlpha", "refStartOther",
       "refDTmin", "refNonPosH", "refLowP", "refDrift", "refAlpha", "refCap",
       "maxSubstepsOneUpdate", "lastSubsteps", "lastRefuseCode", "maxAlphaRatio",
       "lastAlphaRatio", "lastF", "entryOverKappa", "rejectedReversal"]
ME = ["updates", "meCalls", "substeps", "accepted", "rejectedErr", "forcedAtDTmin",
      "forcedClampMc", "rejectedLowP", "abandonedLowP", "capHits", "entryPminClamps",
      "pnResets", "maxSubstepsOneUpdate", "lastSubsteps", "lastForcedAtDTmin",
      "lastAbandonedLowP", "lastCapHit"]
REFCODES = {16: "startOutsideYield", 17: "startAlphaOutsideBounding",
            18: "startInadmissible", 19: "errorAtDTmin", 20: "loadingNonPosH",
            21: "tensionAtDTmin(lowP)", 22: "driftFailed", 23: "alphaOutsideAtDTmin",
            24: "maxSubsteps"}
REFLINE = re.compile(r"refusals step (\d+): (\d+) by code (\{[^}]*\}); n points (\d+); top: (.*)")


def steps(leg):
    return list(csv.DictReader(open(os.path.join(R, leg, "steps.csv"))))


def reflines(leg):
    out = {}
    for ln in open(os.path.join(R, leg, "logs", "log.log"), errors="replace"):
        m = REFLINE.search(ln)
        if m:
            out[int(m.group(1))] = (int(m.group(2)), m.group(3), int(m.group(4)), m.group(5)[:160])
    return out


def fails_by_step(leg):
    recs = rung_fail.parse(leg)
    per = {}
    for r in recs:
        per.setdefault(r["step"], []).append(
            f"{r['rung']}:{r['cause'][0]}" + (f"({float(r['norm_over_tol']):.1f})" if r["norm_over_tol"] else ""))
    return recs, per


summary = {}
os.makedirs(os.path.join(A, "tables"), exist_ok=True)
for leg in LEGS:
    S = steps(leg)
    recs, per = fails_by_step(leg)
    RL = reflines(leg)
    nlast = int(S[-1]["step"])
    # ---- last 20 steps table
    rows = []
    for r in S[-20:]:
        k = int(r["step"])
        rl = RL.get(k)
        rows.append(dict(step=k, s_over_B=f"{float(r['s_over_B']):.6f}", q_kPa=f"{float(r['q_kPa']):.2f}",
                         ds_m=r["ds_m"], rung="NLK"[int(r["rung"])], iters=r["iters"],
                         fails_before=r["fails_before"], failed_rungs=" ".join(per.get(k, [])),
                         refusals=r["cap_step"], ref_codes=(rl[1] if rl else ""),
                         ref_where=(rl[3] if rl else ""), rejErr=r["rejErr_step"],
                         dtmin=r["dtmin_step"], min_p=r["min_p"], max_rho=r["max_rho_alpha"],
                         n_rho_gt_1=r["n_rho_gt_1"], wall_s=r["wall_step_s"]))
    # the failing (floor) step
    rows.append(dict(step=nlast + 1, s_over_B="FLOOR", failed_rungs=" ".join(per.get(nlast + 1, []))))
    with open(os.path.join(A, "tables", f"last20_{leg}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    # ---- rung / cause totals
    n = len(S)
    rung_acc = {c: sum(1 for r in S if r["rung"] == str(i)) for i, c in enumerate("NLK")}
    dsK = sum(float(r["ds_m"]) for r in S if r["rung"] == "2")
    dsT = sum(float(r["ds_m"]) for r in S)
    kits = [int(r["iters"]) for r in S if r["rung"] == "2"]
    last20 = S[-20:]
    causes = {}
    for r in recs:
        causes[f"{r['rung']}:{r['cause']}"] = causes.get(f"{r['rung']}:{r['cause']}", 0) + 1
    floor = [r for r in recs if r["step"] == nlast + 1]
    fc = {}
    for r in floor:
        fc[r["cause"]] = fc.get(r["cause"], 0) + 1
    # N-rung residual ratio at failure in the last 40 steps
    nrat = [float(r["norm_over_tol"]) for r in recs if r["rung"] == "N" and r["norm_over_tol"]
            and r["step"] > nlast - 40]
    # refusals over the leg
    first_K_streak = None
    run = 0
    for r in S:
        run = run + 1 if r["rung"] == "2" else 0
        if run == 10 and first_K_streak is None:
            first_K_streak = (int(r["step"]) - 9, float(r["s_over_B"]))
    summary[leg] = dict(steps=n, s_over_B_end=float(S[-1]["s_over_B"]), q_end=float(S[-1]["q_kPa"]),
                        accepted_by_rung=rung_acc, frac_K=rung_acc["K"] / n,
                        frac_K_last20=sum(1 for r in last20 if r["rung"] == "2") / len(last20),
                        frac_settlement_on_K=dsK / dsT,
                        K_iters_median=(float(np.median(kits)) if kits else None),
                        K_accepted_at_it1=sum(1 for i in kits if i == 1),
                        first_10_K_in_a_row=first_K_streak,
                        failed_rungs_total=len(recs), failed_rung_causes=causes,
                        floor_ladder_rungs=len(floor), floor_ladder_causes=fc,
                        N_fail_norm_over_tol_last40=(dict(n=len(nrat), median=float(np.median(nrat)),
                                                          min=float(min(nrat)), max=float(max(nrat)))
                                                     if nrat else None),
                        refusals_total=sum(int(r["cap_step"]) for r in S),
                        refusals_last20=sum(int(r["cap_step"]) for r in last20),
                        ds_last=float(S[-1]["ds_m"]),
                        ds_min_accepted=min(float(r["ds_m"]) for r in S))

# ---- final-ladder census diff (codes + locations of the refusers at the floor)
LASTCK = {"E_A": 185, "E_B": 375, "E_D": 385, "E_B16": 75}
fl_rows = []
for leg, ck in LASTCK.items():
    z0 = np.load(os.path.join(CK, leg, "ckpt", f"field_step{ck:05d}.npz"))
    z1 = np.load(os.path.join(CK, leg, "ckpt", "field_last_converged.npz"))
    d = z1["stats"] - z0["stats"]
    names = ME if z1["stats"].shape[1] == len(ME) else SAS
    S = steps(leg)
    between = [r for r in S if int(r["step"]) > ck]
    ref_between = sum(int(r["cap_step"]) for r in between)
    info = dict(ckpt=ck, steps_after_ckpt=[int(r["step"]) for r in between],
                refusals_in_those_converged_steps=ref_between)
    if names is SAS:
        tot = int(d[:, 15].sum())
        info["refusals_ckpt_to_floor_incl_floor_ladder"] = tot
        info["floor_ladder_refusals_est"] = tot - ref_between
        info["by_code_ckpt_to_floor"] = {REFCODES[i]: int(d[:, i].sum()) for i in REFCODES if d[:, i].sum() > 0}
        refcol = 15
    else:
        tot = int(d[:, 9].sum())
        info["capHits_ckpt_to_floor_incl_floor_ladder"] = tot
        info["floor_ladder_capHits_est"] = tot - ref_between
        info["forcedAtDTmin_ckpt_to_floor"] = int(d[:, 5].sum())
        info["rejectedLowP_ckpt_to_floor"] = int(d[:, 7].sum())
        info["abandonedLowP_ckpt_to_floor"] = int(d[:, 8].sum())
        refcol = 9
    ks = np.where(d[:, refcol] > 0)[0]
    info["n_points_refusing"] = int(len(ks))
    p, rho, eta = z1["p"], z1["rho"], z1["eta"]
    for k in ks[np.argsort(-d[ks, refcol])]:
        row = dict(leg=leg, k=int(k), element=int(z1["tag"][k]), gp=int(k % 4 + 1),
                   x=round(float(z1["gx"][k]), 3), y=round(float(z1["gy"][k]), 3),
                   n_ref=int(d[k, refcol]), p_kPa=round(float(p[k]), 3),
                   eta=round(float(eta[k]), 3), rho_alpha=round(float(rho[k]), 4))
        if names is SAS:
            row["codes"] = "+".join(f"{REFCODES[i]}x{int(d[k, i])}" for i in REFCODES if d[k, i] > 0)
        else:
            row["codes"] = f"capHits x{int(d[k, 9])}; forcedAtDTmin x{int(d[k, 5])}"
        fl_rows.append(row)
    # where the rho>1 points are
    o = np.where(rho > 1.0)[0]
    info["n_rho_gt_1"] = int(len(o))
    if len(o):
        info["rho_gt_1_x_range"] = [float(z1["gx"][o].min()), float(z1["gx"][o].max())]
        info["rho_gt_1_y_range"] = [float(z1["gy"][o].min()), float(z1["gy"][o].max())]
        info["rho_gt_1_p_range"] = [float(p[o].min()), float(p[o].max())]
        info["rho_max"] = float(np.nanmax(rho))
    info["p_min"] = float(p.min())
    info["n_p_lt_1kPa"] = int((p < 1.0).sum())
    summary[leg]["floor_census"] = info
with open(os.path.join(A, "tables", "floor_refusers.csv"), "w", newline="") as f:
    keys = ["leg", "k", "element", "gp", "x", "y", "n_ref", "codes", "p_kPa", "eta", "rho_alpha"]
    w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(fl_rows)
json.dump(summary, open(os.path.join(A, "tables", "walls_summary.json"), "w"), indent=1)

# ---- step-size history plot
fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
col = dict(E_A="C3", E_B="C0", E_D="C2", E_B16="C1", E_C2="C4")
for leg in LEGS:
    S = steps(leg)
    x = [float(r["s_over_B"]) for r in S]
    ax[0].semilogy(x, [float(r["ds_m"]) for r in S], ".-", ms=3, lw=.7, color=col[leg], label=leg)
    kk = [(float(r["s_over_B"]), int(r["fails_before"])) for r in S if int(r["fails_before"]) > 0]
    ax[1].plot(x, np.cumsum([int(r["cap_step"]) for r in S]), color=col[leg], label=f"{leg} refusals (cum)")
ax[0].axhline(2e-7, color="k", ls=":", label="floor 2e-7 m")
ax[0].axvspan(0.026, 0.041, color="grey", alpha=.15, label="TIMs ME wall")
ax[1].axvspan(0.026, 0.041, color="grey", alpha=.15)
ax[0].set_ylabel("accepted ds (m)"); ax[1].set_ylabel("cumulative refusals / capHits")
ax[1].set_yscale("symlog"); ax[1].set_xlabel("s/B"); ax[0].legend(fontsize=7, ncol=3); ax[1].legend(fontsize=7)
ax[0].set_title("Step-size history and cumulative material refusals")
fig.tight_layout(); fig.savefig(os.path.join(A, "ds_history.png"), dpi=130)

# ---- rung per step (last 60) for the finished legs
fig, axs = plt.subplots(4, 1, figsize=(10, 9))
for a, leg in zip(axs, ["E_A", "E_B", "E_D", "E_B16"]):
    S = steps(leg)[-60:]
    st = [int(r["step"]) for r in S]
    a.bar(st, [int(r["iters"]) for r in S], color=[("C0", "C1", "C3")[int(r["rung"])] for r in S])
    a2 = a.twinx(); a2.semilogy(st, [float(r["ds_m"]) for r in S], "k.-", ms=3, lw=.6)
    a.set_ylabel("iters"); a2.set_ylabel("ds")
    a.set_title(f"{leg}: last 60 steps: bar colour = accepting rung (blue N, orange LS, red Krylov@10x)", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(A, "last60_rungs.png"), dpi=120)
print(json.dumps(summary, indent=1))
