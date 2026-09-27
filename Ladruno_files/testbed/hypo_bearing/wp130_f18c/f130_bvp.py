"""WP-130 / TIMs F18(c): F12's bearing deck with the new CPPM flags, instrumented.

Same idea as `../adr92_f12/f12_bvp.py` (load `sanisand_tau0_band.py` as a module,
force its `INT_SCHEME`, call its own `main()`, so no repo driver is edited and
its controls / JSON provenance still run), plus two seams F12 did not have:

  * `ops.nDMaterial` is wrapped so the LadrunoSANISAND command gets the extra
    WP-130 flags given after `--`, e.g.
        -- -cppmOnFail refuse -cppmHalvings 0 -cppmLineSearch on
  * `ops.analyze` / `ops.algorithm` are wrapped so EVERY analyze call (every
    rung of every attempt) is logged to `<out>/analyze_log.jsonl` with the
    rung's algorithm, the return code, `testIter()`, `testNorm()` (the
    per-iteration unbalance norms -- the convergence-rate evidence), the time
    and the wall seconds it cost.

At the end of the leg the per-Gauss-point `substepStats` census is summed over
the mesh and written to `<out>/census.json`.

usage (from this folder, with the WP-130 build):
    LADRUNO_DIST_BIN=<wt>/dist/bin python -u f130_bvp.py <IntScheme> --out <dir> \
        --legs h1.0_e0.6944 --xlim 10 --zbot 8 --wall 1200 --maxsubsteps 20000 \
        --tantype 2 [-- <extra LadrunoSANISAND flags>]
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

_HYPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DRIVER = os.environ.get("F12_DRIVER", os.path.join(_HYPO, "sanisand_tau0_band.py"))
os.environ.setdefault("LADRUNO_A2_EXPECT_BUILD", "any")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

scheme = int(sys.argv[1])
argv = sys.argv[2:]
extra = []
if "--" in argv:
    k = argv.index("--")
    extra = argv[k + 1:]
    argv = argv[:k]
out_dir = argv[argv.index("--out") + 1]
os.makedirs(out_dir, exist_ok=True)

spec = importlib.util.spec_from_file_location("tau0band", DRIVER)
mod = importlib.util.module_from_spec(spec)
sys.modules["tau0band"] = mod
spec.loader.exec_module(mod)
mod.INT_SCHEME = scheme
ops = mod.ops

_extra_conv = []
for tok in extra:
    try:
        _extra_conv.append(int(tok))
    except ValueError:
        try:
            _extra_conv.append(float(tok))
        except ValueError:
            _extra_conv.append(tok)

_orig_nd = ops.nDMaterial
_orig_an = ops.analyze
_orig_al = ops.algorithm
_state = {"algo": None, "n": 0}
_log = open(os.path.join(out_dir, "analyze_log.jsonl"), "w")


def _nd(*a):
    if len(a) > 1 and a[0] == "LadrunoSANISAND":
        a = tuple(a) + tuple(_extra_conv)
        print("@@F130 nDMaterial", " ".join(str(x) for x in a), flush=True)
    return _orig_nd(*a)


def _al(*a):
    _state["algo"] = a[0] if a else None
    return _orig_al(*a)


def _an(*a):
    t = time.time()
    rc = _orig_an(*a)
    dt = time.time() - t
    try:
        norms = list(ops.testNorm())
    except Exception:
        norms = []
    rec = dict(k=_state["n"], algo=_state["algo"], rc=int(rc), iters=int(ops.testIter()),
               norms=norms, time=ops.getTime(), wall=dt)
    _state["n"] += 1
    _log.write(json.dumps(rec) + "\n")
    _log.flush()
    return rc


ops.nDMaterial = _nd
ops.analyze = _an
ops.algorithm = _al
print("@@F130 INT_SCHEME forced to %d; extra flags %s" % (scheme, _extra_conv), flush=True)

rc = 1
try:
    rc = mod.main(argv)
finally:
    # census over the mesh (the model is still in memory after main returns)
    tot = None
    try:
        for e in ops.getEleTags():
            for gp in range(1, 9):
                try:
                    s = list(ops.eleResponse(e, "material", gp, "substepStats"))
                except Exception:
                    s = []
                if not s:
                    continue
                tot = s if tot is None else [x + y for x, y in zip(tot, s)]
    except Exception as exc:  # pragma: no cover
        print("@@F130 census failed:", exc, flush=True)
    names = ["updates", "meCalls", "substeps", "accepted", "rejectedErr", "forcedAtDTmin",
             "forcedClampMc", "rejectedLowP", "abandonedLowP", "capHits", "entryPminClamps",
             "pnResets", "maxSubstepsOneUpdate", "lastSubsteps", "lastForcedAtDTmin",
             "lastAbandonedLowP", "lastCapHit", "cppmCalls", "cppmNewtonFail",
             "cppmHalvings", "cppmExplicitFail", "cppmExplicitLowP", "cppmRefusals",
             "meFallbacks", "meFallbackOk", "lastCppmRefused"]
    with open(os.path.join(out_dir, "census.json"), "w") as fh:
        json.dump(dict(scheme=scheme, extra=_extra_conv,
                       build=ops.ladrunoBuild() if hasattr(ops, "ladrunoBuild") else "?",
                       census=dict(zip(names, tot)) if tot else None), fh, indent=1)
    _log.close()
raise SystemExit(rc)
