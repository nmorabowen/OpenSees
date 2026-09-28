"""WP-128 Q3: markdown table of out/q3_worst_point.json for the b8 worst point
(1950/3).  Cell = substeps / returned eta / returned f, with flags:
  C  the dT_min forced accept fired with the Mc clamp (eta teleported to Mc)
  F  a forced accept without clamp;  X  rc != 0;  +  f_after > 1e-6 (outside
  the yield surface, reported as success when rc = 0).
CPPM cells: substeps > 0 means CPPM FELL BACK to ModifiedEuler (its own
Newton failed; the fallback runs ME at the prototype's TolR).
Output: out/q3_table.md
"""
import json

import _boot as B

KEY = "(1950, 3)"
INTEG = ["ME(campaign)", "ME(1e-8,honor)", "RK45(1e-10)", "CPPM(1e-10)", "ME+aErr(port)"]


def main():
    with open(f"{B.OUT}/q3_worst_point.json") as fh:
        res = json.load(fh)[KEY]["runs"]
    probes = []
    for k in res:
        n, p = k.split("|")
        if p not in probes:
            probes.append(p)
    lines = ["| probe | " + " | ".join(INTEG) + " |", "|---|" + "---|" * len(INTEG)]
    for p in probes:
        cells = []
        for n in INTEG:
            r = res.get(f"{n}|{p}")
            if r is None:
                cells.append("-")
                continue
            flag = ""
            clamp = r.get("clampMc", 0)
            forced = r.get("forced", 0)
            if clamp or abs(r["eta"] - B.MC) < 1e-9:   # RK45 is not instrumented: clamp signature
                flag += "C"
            elif forced:
                flag += "F"
            if r["rc"] != 0:
                flag += "X"
            if r["f_after"] > 1e-6:
                flag += "+"
            sub = r.get("substeps", 0)
            if n.startswith("RK45"):
                sub = "n/a"
            cells.append(f"{sub} / {r['eta']:.3g} / {r['f_after']:.2g} {flag}".strip())
        lines.append(f"| {p} | " + " | ".join(cells) + " |")
    txt = "\n".join(lines)
    print(txt)
    with open(f"{B.OUT}/q3_table.md", "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
