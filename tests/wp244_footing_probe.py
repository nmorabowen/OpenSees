"""WP-244 runner (fix #1112): a footing case with the footing box's lower bound at -lead * h
(still `track y`), so the old box (lead 0) and the fix (lead 0.5) run on one build. It writes the
records under tests/data/mpm/wp244/ (results: Ladruno_implementation/244_mpm_footing_track.md).

  python3.12 tests/wp244_footing_probe.py s3 --tests tests --mat hs --r 2 --lead 0.5 --s-end 0.5       --out tests/data/mpm/wp244/f1_hs_r2_lead0.5.json                       # F1 (S3 harness)
  python3.12 tests/wp244_footing_probe.py l4 --tests tests --case footing_mildcone_r2_al       --lead 0.0 --steps 210 --out X.json                                    # F3 (L4 ladder)
  ... l4 --case footing_mildcone_r2_al --al-off ...                          # 232f's footing
The engine is the opensees module on PYTHONPATH (dist/bin); the records carry ladrunoBuild.
Records: OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 (results move ~1e-13 with the thread count).
"""
import argparse
import json
import sys
import time
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["s3", "l4"])
    ap.add_argument("--tests", required=True)
    ap.add_argument("--mat", default="hs")
    ap.add_argument("--deck", default=None)
    ap.add_argument("--r", type=int, default=2)
    ap.add_argument("--lead", type=float, default=0.5)
    ap.add_argument("--s-end", type=float, default=0.6)
    ap.add_argument("--case")
    ap.add_argument("--cases-file", default=None, help="tests/mpm_implicit_cases.txt for 232f")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--al-off", action="store_true", help="anti-locking off (232f's footing)")
    ap.add_argument("--wall-cap", type=float, default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.tests)
    from _testbed import ops
    t0 = time.time()
    if a.what == "s3":
        import ladruno_mpm_soil_s_harness as S
        h = 0.5 / 2 ** (a.r - 1)

        class Proxy:                                  # the footing box's lo_y -> -lead*h
            def __getattr__(self, k):
                return getattr(ops, k)

            def mpm(self, *args):
                if len(args) > 4 and args[0] == "gridvel" and args[1] == "footing":
                    args = list(args)
                    args[4] = -a.lead * h
                return ops.mpm(*args)
        pops = Proxy()
        n = 4 * 2 ** (a.r - 1)
        if a.deck:
            ref = S.load(f"s3_fe_{a.mat}_{a.deck}_n{n}.json")
            extra = [S._num(str(x)) for x in ref["san_opts"]]
        else:
            ref = S.load(f"s3_fe_{a.mat}_n{n}.json")
            extra = ()
        res = S.s3_run(pops, ref["decks"][a.mat]["cmd"], a.r, "fbar", extra_opts=extra,
                       s_end=a.s_end, log=lambda m: print(m, flush=True), wall_cap=a.wall_cap)
        res["fe_q_at"] = ref.get("q_at")
    else:
        import ladruno_mpm_antilocking_harness as A
        cs = A.parse_cases(Path(a.cases_file)) if a.cases_file else A.parse_cases()
        c = next(x for x in cs if x["name"] == a.case)
        h = c["h"]
        for g in c["gridvels"]:
            if g[0] == "footing":
                g[2 + 1] = repr(-a.lead * h)        # lo_y token (box lo_x lo_y hi_x hi_y)
        if a.steps:
            c["steps"] = a.steps
        if a.al_off:
            c["antilocking"] = "off"
        r = A.run_ladder(c, ops, log=A.flushprint)
        res = {"antilocking": c.get("antilocking")}
        res.update({k: r.get(k) for k in ("failed_attempts", "relaxed_steps", "halvings", "stopped",
                                     "zone_fallback_share")})
        rr = A.run_record(c, r)
        res.update(t=rr["t"], q=rr["q"], s=rr.get("s"), verdicts=rr.get("verdicts"))
        res["gridvel"] = [g for g in c["gridvels"] if g[0] == "footing"]
    res.update(lead=a.lead, build=str(ops.ladrunoBuild()), wall_s=time.time() - t0)
    Path(a.out).write_text(json.dumps(res, indent=1), encoding="utf-8")
    print("DONE", a.out, res.get("stopped"), res.get("failed_attempts"), flush=True)


if __name__ == "__main__":
    main()
