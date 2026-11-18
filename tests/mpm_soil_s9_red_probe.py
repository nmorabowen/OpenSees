"""WP-232k S9: the demonstrated red of the Zone-A soil oracle (test_s9_onepoint_triax).

Runs the oracle's legs (hs_ps_triax, sanisand_ps_triax; modes A and F) as the test does, then
under three mutants, and reports the first step at which the oracle's rule fails:
  no_flip      -- `mpm stage 1 1` is never issued (the stage flip skipped at particles);
  skip_update  -- the particles' material update is skipped on step 11, the first stage-1 step
                  (that step's grid motion is withheld, so the material sees no increment; the
                  run then continues on the record's increments);
  stale_stress -- the particle stress is not refreshed by the update: each step reports the
                  previous step's committed stress (the update's output dropped).
A mutant is caught when the oracle fails; the unmutated run must pass.

    python tests/mpm_soil_s9_red_probe.py [--out <json>]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ladruno_mpm_soil_s_harness as S  # noqa: E402


class _Ops:
    """ops proxy applying one mutant to the S1 driver's calls."""

    def __init__(self, ops, mutant):
        self._ops, self._m, self._step, self._prev = ops, mutant, 0, None

    def __getattr__(self, name):
        return getattr(self._ops, name)

    def mpm(self, *a):
        if self._m == "no_flip" and a[:1] == ("stage",) and len(a) == 3 and a[2] == 1:
            return None
        if a[:1] == ("gridvel",) and len(a) > 2 and a[2] == "setAffine":
            self._step += 1
            if self._m == "skip_update" and self._step == 11:
                return self._ops.mpm(*a[:3], *([0.0] * (len(a) - 3)))
        if self._m == "stale_stress" and a[:2] == ("particleData", "sigHat"):
            cur = self._ops.mpm(*a)
            out, self._prev = (self._prev if self._prev is not None else cur), cur
            return out
        return self._ops.mpm(*a)


def judge(out):
    bad = [r for r in out["steps"] if not (r["err_sig"] <= r["tol"] and r["err_state"] <= r["tol"])]
    return {"passed": not bad, "first_fail": bad[0] if bad else None, "max_err": out["max_err"]}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    from _testbed import ops
    res = {"build": ops.ladrunoBuild(), "rows": {}}
    for leg in ("hs_ps_triax", "sanisand_ps_triax"):
        for mode in ("A", "F"):
            for mutant in ("none", "no_flip", "skip_update", "stale_stress"):
                o = ops if mutant == "none" else _Ops(ops, mutant)
                j = judge(S.s1_run(o, leg, mode))
                res["rows"][f"{leg}:{mode}:{mutant}"] = j
                ff = j["first_fail"]
                print(f"{leg:18s} {mode} {mutant:12s} "
                      + ("PASS" if j["passed"] else
                         f"RED at step {ff['step']} (err_sig {ff['err_sig']:.3g}, "
                         f"err_state {ff['err_state']:.3g}, tol {ff['tol']:.0e})"))
    caught = all(not v["passed"] for k, v in res["rows"].items() if not k.endswith(":none"))
    clean = all(v["passed"] for k, v in res["rows"].items() if k.endswith(":none"))
    res.update(all_mutants_caught=caught, unmutated_pass=clean)
    print("all mutants caught:", caught, "| unmutated pass:", clean)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1)
    return 0 if (caught and clean) else 1


if __name__ == "__main__":
    sys.exit(main())
