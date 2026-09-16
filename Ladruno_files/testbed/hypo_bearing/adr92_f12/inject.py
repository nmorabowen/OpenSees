"""Inject the generated tables into F12_phase_a_report.md's placeholders."""
import io
import json
import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def section(text, name):
    """Pull one '## X. ...' block out of analyse.py's output."""
    lines = text.splitlines()
    out, on = [], False
    for ln in lines:
        if ln.startswith("## "):
            on = ln.startswith("## " + name)
            continue
        if on:
            out.append(ln)
    return "\n".join(out).strip()


def ctrl_table():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "data", "ct_*.json"))):
        d = json.load(open(f, encoding="utf-8"))
        rows.append((os.path.basename(f)[:-5], d))
    out = ["| run | scheme | N | dEz | global tol | steps done | stalled@ | "
           "ms/step (passed) | wall of the FAILING step (s) | it/step | it max | "
           "subME max | eta_end |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for tag, d in rows:
        out.append("| `%s` | %d | %d | %.1e | %.0e | %d | %s | %.1f | %.2f | %.2f "
                   "| %d | %.0f | %.4f |"
                   % (tag, d["scheme"], d["nstep"], d["ez_max"] / d["nstep"],
                      d.get("push_tol", 0), d["steps_done"], d["stalled_at"],
                      d["ms_per_step_ok"], d["wall_failed_step"], d["iters_mean"],
                      d["iters_max"], d["subME_max"], d["eta_end"]))
    return "\n".join(out)


def main():
    tables = subprocess.run([sys.executable, os.path.join(HERE, "analyse.py")],
                            capture_output=True, text=True).stdout
    rep = os.path.join(HERE, "F12_phase_a_report.md")
    s = open(rep, encoding="utf-8").read()
    s = s.replace("<!--TABLE-A-->",
                  section(tables, "A."))
    s = s.replace("<!--TABLE-B-->", section(tables, "B."))
    s = s.replace("<!--TABLE-D-->", section(tables, "D."))
    s = s.replace("<!--TABLE-E-->", section(tables, "E."))
    s = s.replace("<!--TABLE-C2-->", ctrl_table())
    pb = os.path.join(HERE, "phase_b.md")
    if os.path.exists(pb):
        s = s.replace("<!--PHASE-B-->", open(pb, encoding="utf-8").read().strip())
    open(rep, "w", encoding="utf-8").write(s)
    print("injected;", "PHASE-B" in s and "placeholder still open" or "phase b in")


if __name__ == "__main__":
    main()
