"""ADR-92 P0 / gate G0 — the BINARY half.

Drives a single-Gauss-point drained triaxial path on `LadrunoSANISAND` through the
fork binary and dumps, per step, the full strain (6), stress (6) and the class's
26-entry `state` vector. The numpy oracle replays this EXACT strain sequence; if it
cannot reproduce the stress to 1e-8 the oracle is measuring itself and P0 stops.

Instrument adapted from the TIMs act's own `labs/sfim/probe_sanisand_tx.tcl` (same
cube, same face-load consolidation, same SP-pattern + LoadControl push) so the two
campaigns are reading the same thing. Material constants are the D-L cell's, not
Toyoura's — see README section 2 for provenance.

Run with the installed build's interpreter:

    C:/Users/nmb/venv/opensees_env/Scripts/python.exe probe_binary_triaxial.py --p0 100

Every run opens with `ops.ladrunoBuild()`; the hash goes in the CSV header, because a
number from an unknown build is not evidence.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import opensees as ops

# --- the D-L cell's constants (README section 2; retyped by TIMs from oracle_sheet.yaml) ---
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


def build(p0: float, e_init: float, scheme: int, tan_type: int, tol: float) -> None:
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
    ops.nDMaterial("LadrunoSANISAND", 1, *args,
                   scheme, tan_type, 1, tol, tol,
                   "-Presidual", PRESIDUAL, "-Pmin", PMIN, "-honorTolR", HONOR_TOLR)
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


def consolidate() -> None:
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-8, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    if ops.analyze(10) != 0:
        raise SystemExit("G0: elastic consolidation FAILED")
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    # Re-equilibrate at CONSTANT load. The Tcl instrument this is adapted from keeps
    # LoadControl 0.1 here and takes 5 more steps, so its "p0 = 100" cube actually
    # consolidates to 150 kPa (time reaches 1.5). Held at 0.0 the label is true.
    ops.integrator("LoadControl", 0.0)
    if ops.analyze(5) != 0:
        raise SystemExit("G0: stage-1 re-equilibration FAILED")
    ops.loadConst("-time", 0.0)


def push(nstep: int, ez_max: float, out: str, meta: dict) -> int:
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
            + [f"state{i}" for i in range(26)])
    done = 0
    with open(out, "w", newline="", encoding="utf-8") as fh:
        for k, v in meta.items():
            fh.write(f"# {k}: {v}\n")
        w = csv.writer(fh)
        w.writerow(cols)
        # step 0 = the committed state the push starts from
        for i in range(nstep + 1):
            if i > 0:
                if ops.analyze(1) != 0:
                    print(f"@@stalled at step {i}")
                    break
                done = i
            eps = ops.eleResponse(1, "material", 1, "strain")
            sig = ops.eleResponse(1, "material", 1, "stress")
            st = ops.eleResponse(1, "material", 1, "state")
            if len(eps) != 6 or len(sig) != 6:
                raise SystemExit(f"G0: unexpected response widths {len(eps)}/{len(sig)}")
            st = list(st) + [float("nan")] * (26 - len(st))
            w.writerow([i, ops.getTime()] + list(eps) + list(sig) + st[:26])
    return done


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p0", type=float, default=100.0)
    ap.add_argument("--e-init", type=float, default=CONSTS["e_init"])
    ap.add_argument("--scheme", type=int, default=1, help="1=ModifiedEuler, 2=BackwardEuler")
    ap.add_argument("--tan-type", type=int, default=2, help="2=consistent (see #792 T2)")
    ap.add_argument("--tol", type=float, default=1.0e-10)
    ap.add_argument("--nstep", type=int, default=400)
    ap.add_argument("--ez-max", type=float, default=0.20)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    build_hash = ops.ladrunoBuild()
    out = a.out or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "data",
        f"tx_p{a.p0:g}_e{a.e_init:g}_s{a.scheme}_n{a.nstep}.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    meta = dict(build=build_hash, pyd=ops.__file__, python=sys.version.split()[0],
                p0=a.p0, e_init=a.e_init, scheme=a.scheme, tan_type=a.tan_type,
                tol=a.tol, nstep=a.nstep, ez_max=a.ez_max,
                Pmin=PMIN, Presidual=PRESIDUAL, honorTolR=HONOR_TOLR,
                consts={k: (a.e_init if k == "e_init" else CONSTS[k]) for k in ORDER})

    build(a.p0, a.e_init, a.scheme, a.tan_type, a.tol)
    consolidate()
    done = push(a.nstep, a.ez_max, out, meta)

    print(f"@@build: {build_hash}")
    print(f"@@out: {out}")
    print(f"@@steps: {done}/{a.nstep}")


if __name__ == "__main__":
    main()
