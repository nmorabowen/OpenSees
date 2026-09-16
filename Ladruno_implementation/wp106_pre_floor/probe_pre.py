"""WP-106 / ADR-93 II.1 -- the BINARY probe for the elastic-only floor `-pRe`.

A single-element drained triaxial on `LadrunoSANISAND`, driven through the fork
binary, dumping per step: strain(6), stress(6), the 26-entry `state` vector, the
`psi` and `yieldDistance` diagnostics and the ADR-86b `substeps` pair.

Adapted from `adr92_p0_oracle/probe_binary_triaxial.py` (same cube, same
face-load consolidation, same SP + LoadControl push) so the three campaigns read
the same instrument; the only addition is `--pRe`.

    python3.12 probe_pre.py --p0 20 --pRe 0   --out data/tx_p20_pre0.csv
    python3.12 probe_pre.py --p0 20 --pRe 1.0 --out data/tx_p20_pre1.csv

`--bin` points `sys.path` at a `dist/bin` so the SAME script can drive the
pre-change binary (which has no `-pRe` flag -- ask for pRe 0 there).
Every run opens with `ops.ladrunoBuild()`: a number from an unknown build is not
evidence.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys

ap0 = argparse.ArgumentParser(add_help=False)
ap0.add_argument("--bin", default=None)
_known, _rest = ap0.parse_known_args()
if _known.bin:
    sys.path.insert(0, os.path.abspath(_known.bin))

import opensees as ops  # noqa: E402

# the D-L cell's constants (adr92_p0_oracle/README section 2)
CONSTS = dict(
    G0=264.32, nu=0.3129, e_init=0.6944, Mc=1.3309, c=0.71, lambda_c=0.027,
    e0=0.83, ksi=0.45, P_atm=101.0, m=0.005, h0=1.3, ch=0.968, nb=3.5,
    A0=0.05, nd=5.75, z_max=12.5, cz=1100.0, rho=2.0,
)
ORDER = ["G0", "nu", "e_init", "Mc", "c", "lambda_c", "e0", "ksi", "P_atm", "m",
         "h0", "ch", "nb", "A0", "nd", "z_max", "cz", "rho"]

PMIN = 0.0101
PRESIDUAL = 0.0
HONOR_TOLR = 0


def build(p0, e_init, scheme, tan_type, tol, pre, max_substeps=0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in enumerate(
        [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)], start=1):
        ops.node(tag, float(x), float(y), float(z))
    for n in (1, 4, 5, 8):
        ops.fix(n, 1, 0, 0)
    for n in (1, 2, 5, 6):
        ops.fix(n, 0, 1, 0)
    for n in (1, 2, 3, 4):
        ops.fix(n, 0, 0, 1)

    vals = dict(CONSTS)
    vals["e_init"] = e_init
    args = [vals[k] for k in ORDER]
    flags = ["-Presidual", PRESIDUAL, "-Pmin", PMIN, "-honorTolR", HONOR_TOLR]
    if max_substeps:
        flags += ["-maxSubsteps", int(max_substeps)]
    if pre != 0.0:
        flags += ["-pRe", float(pre)]
    ops.nDMaterial("LadrunoSANISAND", 1, *args, scheme, tan_type, 1, tol, tol, *flags)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1, "-formulation", "bbar")

    q4 = -p0 / 4.0
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, q4, 0.0, 0.0)
    for n in (3, 4, 7, 8):
        ops.load(n, 0.0, q4, 0.0)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.0, 0.0, q4)


def consolidate():
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-8, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    if ops.analyze(10) != 0:
        raise SystemExit("probe: elastic consolidation FAILED")
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.integrator("LoadControl", 0.0)
    if ops.analyze(5) != 0:
        raise SystemExit("probe: stage-1 re-equilibration FAILED")
    ops.loadConst("-time", 0.0)


def _resp(name):
    try:
        v = ops.eleResponse(1, "material", 1, name)
    except Exception:
        return []
    return list(v) if v else []


def push(nstep, ez_max, out, meta):
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for n in (5, 6, 7, 8):
        ops.sp(n, 3, -1.0)

    du = ez_max / nstep
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-9, 200, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", du)
    ops.analysis("Static")

    cols = (["step", "pseudo_time"]
            + [f"eps{i}" for i in range(6)]
            + [f"sig{i}" for i in range(6)]
            + [f"state{i}" for i in range(26)]
            + ["psi", "yieldDistance", "substeps_me", "substeps_capHit"])
    done = 0
    with open(out, "w", newline="", encoding="utf-8") as fh:
        for k, v in meta.items():
            fh.write(f"# {k}: {v}\n")
        w = csv.writer(fh)
        w.writerow(cols)
        for i in range(nstep + 1):
            if i > 0:
                if ops.analyze(1) != 0:
                    print(f"@@stalled at step {i}")
                    break
                done = i
            eps = _resp("strain")
            sig = _resp("stress")
            st = list(_resp("state"))
            st = st + [float("nan")] * (26 - len(st))
            psi = _resp("psi")
            yd = _resp("yieldDistance")
            sub = _resp("substeps")
            sub = list(sub) + [float("nan")] * (2 - len(sub))
            w.writerow([i, ops.getTime()] + list(eps) + list(sig) + st[:26]
                       + [psi[0] if psi else float("nan"),
                          yd[0] if yd else float("nan"),
                          sub[0], sub[1]])
    return done


def main():
    ap = argparse.ArgumentParser(parents=[ap0])
    ap.add_argument("--p0", type=float, default=100.0)
    ap.add_argument("--e-init", type=float, default=CONSTS["e_init"])
    ap.add_argument("--scheme", type=int, default=1)
    ap.add_argument("--tan-type", type=int, default=2)
    ap.add_argument("--tol", type=float, default=1.0e-10)
    ap.add_argument("--nstep", type=int, default=400)
    ap.add_argument("--ez-max", type=float, default=0.20)
    ap.add_argument("--pRe", type=float, default=0.0)
    ap.add_argument("--max-substeps", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    build_hash = ops.ladrunoBuild()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    meta = dict(build=build_hash, pyd=ops.__file__, python=sys.version.split()[0],
                p0=a.p0, e_init=a.e_init, scheme=a.scheme, tan_type=a.tan_type,
                tol=a.tol, nstep=a.nstep, ez_max=a.ez_max, pRe=a.pRe,
                Pmin=PMIN, Presidual=PRESIDUAL, honorTolR=HONOR_TOLR,
                max_substeps=a.max_substeps,
                consts={k: (a.e_init if k == "e_init" else CONSTS[k]) for k in ORDER})

    build(a.p0, a.e_init, a.scheme, a.tan_type, a.tol, a.pRe,
          max_substeps=a.max_substeps)
    consolidate()
    done = push(a.nstep, a.ez_max, a.out, meta)

    print(f"@@build: {build_hash}")
    print(f"@@out: {a.out}")
    print(f"@@steps: {done}/{a.nstep}")


if __name__ == "__main__":
    main()
