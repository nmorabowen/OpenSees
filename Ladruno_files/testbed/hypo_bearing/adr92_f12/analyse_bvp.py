"""F12 phase (b)/(c) rollup: the bearing deck, three arms, matched s/B."""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ARMS = [("s1_T0", "IntScheme 1, TanType 0"),
        ("s1_T2", "IntScheme 1, TanType 2"),
        ("s2_T2", "IntScheme 2, TanType 2")]
LEG = "a2_h1.0_e0.6944"
RUNG_CAPS = (25, 40, 60)          # adr92_bvp_gate.py's own decomposition
CONVERGED_ITERS = 5


def load(arm):
    d = os.path.join(HERE, "bvp", arm)
    js = json.load(open(os.path.join(d, LEG + ".json"), encoding="utf-8"))
    rows = []
    with open(os.path.join(d, LEG + "_curve.csv"), newline="", encoding="utf-8") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            rows = list(csv.DictReader([ln] + list(fh)))
            break
    return js, rows


def q_at(rows, sb):
    prev = None
    for r in rows:
        x = float(r["s_over_B"])
        if x >= sb:
            if prev is None:
                return None
            x0, q0 = float(prev["s_over_B"]), float(prev["q_foot_kPa"])
            f = (sb - x0) / (x - x0)
            return q0 + f * (float(r["q_foot_kPa"]) - q0)
        prev = r
    return None


def w_at(rows, sb):
    prev = None
    for r in rows:
        x = float(r["s_over_B"])
        if x >= sb:
            return float(r["wall_s"])
        prev = r
    return None


def main():
    data = {}
    for arm, _ in ARMS:
        try:
            data[arm] = load(arm)
        except (OSError, ValueError):
            pass

    print("### 7.1 Leg outcome (x10z8, `h1.0_e0.6944`, 624 DOF, 84 hexes, "
          "1200 s budget each)\n")
    print("| arm | scheme | TanType | steps | mode | s/B reached | q at end "
          "(kPa) | wall s | s per step | nfail | nsub | nrelax | CLAMPING | "
          "OutsideBounding |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for arm, label in ARMS:
        if arm not in data:
            print("| %s | -- | -- | RUN MISSING | | | | | | | | | | |" % label)
            continue
        j, _ = data[arm]
        print("| %s | %d | %d | %d | %s | %.5f | %.2f | %.0f | %.2f | %d | %d | %d "
              "| %d | %d |"
              % (label, j["int_scheme"], j["tan_type"], j["steps"], j["mode"],
                 j["s_peak_over_B"], j["q_u"], j["wall_s"],
                 j["wall_s"] / max(j["steps"], 1), j["nfail"], j["nsub"],
                 j["nrelax"], j["n_clamping"], j["n_outside_bounding"]))
    print()

    print("### 7.2 Committed load-settlement at matched s/B\n")
    sbs = [0.001, 0.002, 0.005, 0.01, 0.015, 0.02]
    hdr = "| s/B | " + " | ".join(l for _, l in ARMS if _ in data) + " |"
    print(hdr)
    print("|---" * (1 + len([a for a, _ in ARMS if a in data])) + "|")
    for sb in sbs:
        cells = []
        for arm, _ in ARMS:
            if arm not in data:
                continue
            q = q_at(data[arm][1], sb)
            cells.append("%.2f" % q if q is not None else "--")
        print("| %.3f | %s |" % (sb, " | ".join(cells)))
    print()
    if "s1_T2" in data and "s2_T2" in data:
        print("**scheme 2 vs scheme 1 (both TanType 2) at matched s/B — the 1 % bar:**\n")
        print("| s/B | q s1_T2 | q s2_T2 | rel. diff | inside 1 %? |")
        print("|---|---|---|---|---|")
        for sb in sbs:
            a = q_at(data["s1_T2"][1], sb)
            b = q_at(data["s2_T2"][1], sb)
            if a is None or b is None:
                print("| %.3f | %s | %s | -- | (not reached by both) |"
                      % (sb, "%.2f" % a if a else "--", "%.2f" % b if b else "--"))
                continue
            d = abs(b - a) / a
            print("| %.3f | %.2f | %.2f | %.3f %% | %s |"
                  % (sb, a, b, 100 * d, "YES" if d <= 0.01 else "**NO**"))
        print()

    print("### 7.3 Wall clock to matched s/B\n")
    print("| s/B | " + " | ".join(l for a, l in ARMS if a in data) + " |")
    print("|---" * (1 + len([a for a, _ in ARMS if a in data])) + "|")
    for sb in sbs:
        cells = []
        for arm, _ in ARMS:
            if arm not in data:
                continue
            w = w_at(data[arm][1], sb)
            cells.append("%.0f s" % w if w is not None else "not reached")
        print("| %.3f | %s |" % (sb, " | ".join(cells)))
    print()

    print("### 7.4 Ladder decomposition (`adr92_bvp_gate.py`'s own arithmetic)\n")
    print("| arm | steps | rung1 | rung2 | rung3 | past rung 1 % | "
          "failed-rung iteration share % |")
    print("|---|---|---|---|---|---|---|")
    for arm, label in ARMS:
        if arm not in data:
            continue
        j, _ = data[arm]
        st, nf, ns, nr = j["steps"], j["nfail"], j["nsub"], j["nrelax"]
        rung3 = nr
        rung2 = nf - 3 * ns - 2 * nr
        rung1 = st - rung2 - rung3
        past1 = 100.0 * (rung2 + rung3) / max(st, 1)
        fail_it = (ns * sum(RUNG_CAPS) + rung2 * RUNG_CAPS[0]
                   + rung3 * (RUNG_CAPS[0] + RUNG_CAPS[1]))
        ok_it = st * CONVERGED_ITERS
        share = 100.0 * fail_it / max(fail_it + ok_it, 1)
        print("| %s | %d | %d | %d | %d | %.1f | %.1f |"
              % (label, st, rung1, rung2, rung3, past1, share))


if __name__ == "__main__":
    main()
