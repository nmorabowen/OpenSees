"""WP-105 / F12 phase (a) rollup: scheme 1 vs scheme 2 on the replayed strain
path, in the ADR-92 certificate's own metric (relative stress-norm deviation).
"""
import csv
import json
import glob
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")


def load(tag):
    with open(os.path.join(D, tag + ".csv"), newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    with open(os.path.join(D, tag + ".json"), encoding="utf-8") as fh:
        js = json.load(fh)
    return rows, js


def sig(r):
    return [float(r["sig%d" % k]) for k in range(6)]


def nrm(v):
    return math.sqrt(sum(x * x for x in v))


def reldev(a, b):
    """||a-b|| / ||b||, the P0/G0 'sigma rel' metric."""
    return nrm([x - y for x, y in zip(a, b)]) / max(nrm(b), 1e-300)


def resample(rows, n):
    """value of the row at fractional position u in [0,1], by index (paths are
    uniform in pseudo-time)."""
    m = len(rows) - 1
    out = []
    for i in range(n + 1):
        u = i / float(n)
        x = u * m
        j = min(int(round(x)), m)
        out.append(rows[j])
    return out


def main():
    print("## A. Replay certificate -- scheme 2 vs scheme 1 on the SAME strain path\n")
    print("| p0 | N | dEz | scheme | steps | ms/step | it/step | subME>0 | eta_end | "
          "p_end | q_end | rel.dev vs own-scheme fine | rel.dev s2-vs-s1 (same N) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    summary = {}
    for P in (100, 20):
        fine = {}
        for S in (1, 2):
            try:
                fine[S] = load("rp_p%d_s%d_n1820" % (P, S))
            except FileNotFoundError:
                fine[S] = None
        for N in (1820, 182, 40):
            arms = {}
            for S in (1, 2):
                tag = "rp_p%d_s%d_n%d" % (P, S, N)
                try:
                    arms[S] = load(tag)
                except FileNotFoundError:
                    continue
            for S, (rows, js) in sorted(arms.items()):
                own = ""
                if fine.get(S) and N != 1820:
                    fr = resample(fine[S][0], N)
                    n = min(len(rows), len(fr)) - 1
                    own = "%.3e" % max(reldev(sig(rows[i]), sig(fr[i]))
                                       for i in range(1, n + 1))
                cross = ""
                if S == 2 and 1 in arms:
                    r1 = arms[1][0]
                    n = min(len(rows), len(r1)) - 1
                    cross = "%.3e" % max(reldev(sig(rows[i]), sig(r1[i]))
                                         for i in range(1, n + 1))
                dez = 0.0182 / N
                print("| %d | %d | %.1e | %d | %d/%d | %.1f | %.2f | %d | %.4f | %.2f | "
                      "%.2f | %s | %s |"
                      % (P, N, dez, S, js["steps_done"], N, js["ms_per_step_ok"],
                         js["iters_mean"], js["subME_nonzero_steps"], js["eta_end"],
                         js["p_end"], js["q_end"], own, cross))
                summary[(P, N, S)] = js
    print()

    print("## B. terminal-step cross-scheme deviation and eta/M^b\n")
    print("| p0 | N | eta_end s1 | eta_end s2 | d(eta)/eta | M^b (s1) | eta/M^b s1 | "
          "eta/M^b s2 | rel.dev terminal |")
    print("|---|---|---|---|---|---|---|---|---|")
    for P in (100, 20):
        for N in (1820, 182, 40):
            try:
                r1, j1 = load("rp_p%d_s1_n%d" % (P, N))
                r2, j2 = load("rp_p%d_s2_n%d" % (P, N))
            except FileNotFoundError:
                continue
            e1, e2 = j1["eta_end"], j2["eta_end"]
            print("| %d | %d | %.5f | %.5f | %.3e | %.4f | %.4f | %.4f | %.3e |"
                  % (P, N, e1, e2, abs(e2 - e1) / e1, j1["Mb_end"],
                     e1 / j1["Mb_end"], e2 / j2["Mb_end"],
                     reldev(sig(r2[-1]), sig(r1[-1]))))
    print()

    print("## C. free-standing drained-triaxial (global Newton ON) -- robustness\n")
    print("| tag | scheme | N | dEz | done | stalled@ | ms/step(ok) | wall of the "
          "failed step (s) | it/step | it max | subME max | eta peak |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for f in sorted(glob.glob(os.path.join(D, "tx_*.json"))):
        js = json.load(open(f, encoding="utf-8"))
        tag = os.path.basename(f)[:-5]
        mspo = js.get("ms_per_step_ok", js["ms_per_step"])
        wfs = js.get("wall_failed_step", 0.0)
        print("| %s | %d | %d | %.1e | %d | %s | %.1f | %.2f | %.2f | %d | %.0f | %.4f |"
              % (tag, js["scheme"], js["nstep"], js["ez_max"] / js["nstep"],
                 js["steps_done"], js["stalled_at"], mspo, wfs, js["iters_mean"],
                 js["iters_max"], js["subME_max"], js["eta_peak"] or 0.0))
    print()

    print("## D. floor path (p -> p_min, ADR-93 ring regime)\n")
    print("| tag | scheme | N | done | ms/step | subME max | steps with subME>0 | "
          "eta peak | p min | p end |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for f in sorted(glob.glob(os.path.join(D, "fl_*.json"))):
        js = json.load(open(f, encoding="utf-8"))
        tag = os.path.basename(f)[:-5]
        print("| %s | %d | %d | %d | %.2f | %.0f | %d | %.4f | %.4f | %.4f |"
              % (tag, js["scheme"], js["nstep"], js["steps_done"],
                 js.get("ms_per_step_ok", js["ms_per_step"]), js["subME_max"],
                 js["subME_nonzero_steps"], js["eta_peak"] or 0.0,
                 js["p_min_along"], js["p_end"]))
    print()

    print("## E. floor path, step-by-step cross-scheme (n=160)\n")
    try:
        r1, _ = load("fl_p5_s1_n160")
        r2, _ = load("fl_p5_s2_n160")
        print("| step | t | p s1 | p s2 | eta s1 | eta s2 | rel.dev | subME s1 | subME s2 |")
        print("|---|---|---|---|---|---|---|---|---|")
        for i in range(0, min(len(r1), len(r2)), 8):
            a, b = r1[i], r2[i]
            print("| %s | %.4f | %.5g | %.5g | %.4f | %.4f | %.3e | %s | %s |"
                  % (a["step"], float(a["t"]), float(a["p"]), float(b["p"]),
                     float(a["eta"]), float(b["eta"]),
                     reldev(sig(b), sig(a)) if i else 0.0,
                     a["subME"], b["subME"]))
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()
