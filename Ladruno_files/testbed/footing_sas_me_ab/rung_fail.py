"""Classify every FAILED ladder rung in a footing_ab log as DIVERGENCE (the
convergence test ran out of iterations: '<Algo>::solveCurrentStep() -the
ConvergenceTest object failed in test()') or REFUSAL (the material refused an
update: '... -the Integrator failed in update()', incl. the line search's own
'InitialInterpolatedLineSearch::search() -the Integrator failed in update()').
Each failed rung ends with 'StaticAnalysis::analyze() - the Algorithm failed'.
Rungs are N, L, K in order within one attempt (ds); a step line closes the
attempt successfully, 'ladder failed, ds ->' closes it as cut.

NOTE: the SANISAND 'REFUSED (...)' warning prints ONCE PER INTEGRATION POINT
for the whole run, so its code is recorded when present but its absence does
not mean no refusal; the refusal count per step is steps.csv cap_step.

    python rung_fail.py LEG [LEG ...]   -> prints a table, writes runs/<leg>/rung_fail.csv
"""
import csv, os, re, sys
H = os.path.dirname(os.path.abspath(__file__))
STEP = re.compile(r"\] step +(\d+) s/B ([0-9.]+) .* rung (\w) it +(\d+)")
CUT = re.compile(r"step (\d+): ladder failed, ds -> ([0-9.e+-]+)")
CODE = re.compile(r"REFUSED \((\w+), code (\d+)\)")
NORM = re.compile(r"current Norm: ([0-9.e+-]+) \(max: ([0-9.e+-]+)")


def parse(leg):
    lines = open(os.path.join(H, "runs", leg, "logs", "log.log"), errors="replace").read().splitlines()
    recs, cause, codes, rung, attempt, ratio = [], None, set(), 0, 1, None
    nxt = 1
    for ln in lines:
        m = CODE.search(ln)
        if m:
            codes.add(m.group(1))
        m = NORM.search(ln)
        if m:
            ratio = float(m.group(1)) / float(m.group(2))
        if "ConvergenceTest object failed in test()" in ln:
            cause = cause or "DIVERGENCE"
        elif "Integrator failed in update()" in ln:
            cause = "REFUSAL"          # refusal wins over a later test failure
        elif "StaticAnalysis::analyze() - the Algorithm failed" in ln:
            recs.append(dict(step=nxt, attempt=attempt, rung="NLK"[rung % 3],
                             cause=cause or "OTHER", codes="+".join(sorted(codes)),
                             norm_over_tol=(f"{ratio:.3g}" if ratio is not None and cause == "DIVERGENCE" else "")))
            cause, codes, rung, ratio = None, set(), rung + 1, None
        else:
            m = CUT.search(ln)
            if m:
                attempt += 1; rung = 0
                continue
            m = STEP.search(ln)
            if m:
                nxt = int(m.group(1)) + 1; attempt = 1; rung = 0; cause = None; codes = set()
    return recs


def main():
    for leg in sys.argv[1:]:
        recs = parse(leg)
        with open(os.path.join(H, "runs", leg, "rung_fail.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, ["step", "attempt", "rung", "cause", "codes", "norm_over_tol"])
            w.writeheader(); w.writerows(recs)
        tot = {}
        for r in recs:
            tot[(r["rung"], r["cause"])] = tot.get((r["rung"], r["cause"]), 0) + 1
        print(f"== {leg}: failed rungs {len(recs)}; " + ", ".join(
            f"{k[0]}:{k[1]} {v}" for k, v in sorted(tot.items())))
        # per step compact: step -> list of rung:cause
        per = {}
        for r in recs:
            per.setdefault(r["step"], []).append(f"a{r['attempt']}{r['rung']}:{r['cause'][0]}" + (f"({r['norm_over_tol']})" if r['norm_over_tol'] else ""))
        for s in sorted(per)[-int(os.environ.get("NLAST", "12")):]:
            print(f"   step {s}: {' '.join(per[s])}")


if __name__ == "__main__":
    main()
